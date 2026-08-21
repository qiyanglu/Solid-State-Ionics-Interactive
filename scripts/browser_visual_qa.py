"""Capture and inspect local marimo apps through Chromium DevTools.

This lightweight helper is intentionally outside the notebooks. It supports
release-time browser screenshots and simple overlap/clipping warnings without
adding a runtime dependency to the WASM teaching apps.
"""

from __future__ import annotations

import argparse
import base64
import json
import time
from pathlib import Path
from urllib.request import urlopen

import websocket


class CDPClient:
    """Minimal synchronous Chrome DevTools Protocol client."""

    def __init__(self, endpoint: str) -> None:
        self._socket = websocket.create_connection(endpoint, timeout=60)
        self._message_id = 0

    def call(self, method: str, **params):
        self._message_id += 1
        request_id = self._message_id
        self._socket.send(
            json.dumps({"id": request_id, "method": method, "params": params})
        )
        while True:
            message = json.loads(self._socket.recv())
            if message.get("id") == request_id:
                if "error" in message:
                    raise RuntimeError(message["error"])
                return message.get("result", {})

    def evaluate(self, expression: str):
        result = self.call(
            "Runtime.evaluate",
            expression=expression,
            awaitPromise=True,
            returnByValue=True,
        )
        if "exceptionDetails" in result:
            details = result["exceptionDetails"]
            raise RuntimeError(details.get("text", "JavaScript evaluation failed"))
        return result["result"].get("value")

    def close(self) -> None:
        self._socket.close()


def page_endpoint(debug_port: int, target_url: str) -> str:
    with urlopen(f"http://127.0.0.1:{debug_port}/json", timeout=5) as response:
        targets = json.load(response)
    for target in targets:
        if target.get("type") == "page" and target.get("url", "").startswith(target_url):
            return target["webSocketDebuggerUrl"]
    raise RuntimeError(f"No browser tab starts with {target_url!r}")


def wait_until_ready(client: CDPClient, timeout_s: float = 90.0) -> None:
    deadline = time.monotonic() + timeout_s
    while time.monotonic() < deadline:
        ready = client.evaluate(
            "document.readyState === 'complete' && "
            "document.body && document.body.innerText.trim().length > 100"
        )
        if ready:
            time.sleep(2.0)
            return
        time.sleep(0.5)
    raise TimeoutError("The notebook did not finish rendering")


def set_viewport(client: CDPClient, width: int, height: int) -> None:
    client.call(
        "Emulation.setDeviceMetricsOverride",
        width=width,
        height=height,
        deviceScaleFactor=1,
        mobile=width < 600,
    )


def set_range_by_label(client: CDPClient, label_text: str, value: float) -> None:
    script = f"""
    (() => {{
      const labelText = {json.dumps(label_text)};
      const host = [...document.querySelectorAll('marimo-slider')].find((node) =>
        node.getAttribute('data-label')?.includes(labelText));
      if (!host) throw new Error(`No slider labelled ${{labelText}}`);
      const thumb = host.shadowRoot?.querySelector('[role=slider]');
      if (!thumb) throw new Error(`Slider ${{labelText}} has no thumb`);
      const root = thumb.parentElement.parentElement;
      const minimum = Number(thumb.getAttribute('aria-valuemin'));
      const maximum = Number(thumb.getAttribute('aria-valuemax'));
      const clamped = Math.max(minimum, Math.min(maximum, {float(value)}));
      const fraction = (clamped - minimum) / (maximum - minimum);
      const rect = root.getBoundingClientRect();
      return {{x: rect.left + fraction * rect.width, y: rect.top + rect.height / 2}};
    }})()
    """
    point = client.evaluate(script)
    client.call(
        "Input.dispatchMouseEvent",
        type="mousePressed",
        x=point["x"],
        y=point["y"],
        button="left",
        clickCount=1,
    )
    client.call(
        "Input.dispatchMouseEvent",
        type="mouseReleased",
        x=point["x"],
        y=point["y"],
        button="left",
        clickCount=1,
    )
    time.sleep(2.0)


def figure_signature(client: CDPClient) -> str:
    return client.evaluate(
        "Array.from(document.querySelectorAll('img')).map((e) => e.currentSrc).join('|')"
    )


def wait_for_figure_update(
    client: CDPClient, previous_signature: str, timeout_s: float = 120.0
) -> None:
    deadline = time.monotonic() + timeout_s
    last_signature = previous_signature
    stable_since = time.monotonic()
    changed = False
    while time.monotonic() < deadline:
        current_signature = figure_signature(client)
        changed = changed or current_signature != previous_signature
        if current_signature != last_signature:
            last_signature = current_signature
            stable_since = time.monotonic()
        elif changed and time.monotonic() - stable_since >= 8.0:
            return
        time.sleep(1.0)
    raise TimeoutError("Figures did not settle after the control change")


def scroll_to_text(client: CDPClient, text: str) -> None:
    script = f"""
    (() => {{
      const targetText = {json.dumps(text)};
      const candidates = [...document.querySelectorAll('h1, h2, h3, p, figcaption, .marimo-cell')];
      const target = candidates.find((node) =>
        node.textContent && node.textContent.includes(targetText));
      if (!target) throw new Error(`No visible section contains ${{targetText}}`);
      target.scrollIntoView({{block: 'start'}});
      window.scrollBy(0, -24);
      return target.textContent.trim().slice(0, 100);
    }})()
    """
    client.evaluate(script)
    time.sleep(0.7)


def layout_warnings(client: CDPClient):
    script = r"""
    (() => {
      const viewport = {width: innerWidth, height: innerHeight};
      const visible = (element) => {
        const style = getComputedStyle(element);
        const rect = element.getBoundingClientRect();
        return style.display !== 'none' && style.visibility !== 'hidden' &&
          Number(style.opacity) !== 0 && rect.width > 2 && rect.height > 2;
      };
      const figures = [...document.querySelectorAll('img, canvas, svg')].filter(visible);
      const clipped = figures.flatMap((element, index) => {
        const rect = element.getBoundingClientRect();
        const parent = element.parentElement?.getBoundingClientRect();
        const warnings = [];
        if (rect.right > viewport.width + 1 || rect.left < -1) {
          warnings.push(`figure ${index + 1} exceeds viewport horizontally`);
        }
        if (parent && (rect.right > parent.right + 2 || rect.left < parent.left - 2)) {
          warnings.push(`figure ${index + 1} exceeds its container`);
        }
        return warnings;
      });
      const textBlocks = [...document.querySelectorAll('p, h1, h2, h3, label')]
        .filter(visible)
        .map((element) => ({element, rect: element.getBoundingClientRect()}));
      const controls = [...document.querySelectorAll('input, button, select')].filter(visible);
      const overlapped = [];
      for (const control of controls) {
        const a = control.getBoundingClientRect();
        for (const {element, rect: b} of textBlocks) {
          if (control.closest('label')?.contains(element) || element.contains(control)) continue;
          const area = Math.max(0, Math.min(a.right, b.right) - Math.max(a.left, b.left)) *
            Math.max(0, Math.min(a.bottom, b.bottom) - Math.max(a.top, b.top));
          if (area > 12) {
            overlapped.push(`control overlaps text: ${element.textContent.trim().slice(0, 50)}`);
            break;
          }
        }
      }
      return {viewport, figures: figures.length, clipped, overlapped};
    })()
    """
    return client.evaluate(script)


def capture(client: CDPClient, output: Path, full_page: bool) -> None:
    if full_page:
        metrics = client.call("Page.getLayoutMetrics")
        size = metrics["cssContentSize"]
        clip = {
            "x": 0,
            "y": 0,
            "width": size["width"],
            "height": size["height"],
            "scale": 1,
        }
        result = client.call(
            "Page.captureScreenshot",
            format="jpeg",
            quality=85,
            clip=clip,
            captureBeyondViewport=True,
            fromSurface=True,
        )
    else:
        result = client.call(
            "Page.captureScreenshot", format="jpeg", quality=85, captureBeyondViewport=False, fromSurface=True
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(base64.b64decode(result["data"]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--url", required=True)
    parser.add_argument("--port", type=int, default=9222)
    parser.add_argument("--width", type=int, required=True)
    parser.add_argument("--height", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--slider", action="append", default=[], metavar="LABEL=VALUE")
    parser.add_argument("--focus-text")
    parser.add_argument("--scroll-y", type=float)
    parser.add_argument("--scroll-offset", type=float, default=0.0)
    parser.add_argument("--full-page", action="store_true")
    args = parser.parse_args()

    client = CDPClient(page_endpoint(args.port, args.url))
    try:
        client.call("Page.enable")
        client.call("Runtime.enable")
        set_viewport(client, args.width, args.height)
        wait_until_ready(client)
        previous_signature = figure_signature(client) if args.slider else ""
        if args.slider:
            client.evaluate("document.getElementById('App').scrollTop = 0")
            time.sleep(0.5)
        for assignment in args.slider:
            label, raw_value = assignment.rsplit("=", 1)
            set_range_by_label(client, label, float(raw_value))
        if args.slider:
            wait_for_figure_update(client, previous_signature)
        if args.scroll_y is not None:
            client.evaluate(f"document.getElementById('App').scrollTop = {args.scroll_y}")
        elif args.focus_text:
            scroll_to_text(client, args.focus_text)
        else:
            client.evaluate("window.scrollTo(0, 0)")
        if args.scroll_offset:
            client.evaluate(f"document.getElementById('App').scrollTop += {args.scroll_offset}")
        time.sleep(0.5)
        warnings = layout_warnings(client)
        capture_client = CDPClient(page_endpoint(args.port, args.url))
        try:
            capture(capture_client, args.output, args.full_page)
        finally:
            capture_client.close()
        print(
            json.dumps(
                {"screenshot": str(args.output), "layout": warnings}, indent=2
            )
        )
    finally:
        client.close()


if __name__ == "__main__":
    main()
