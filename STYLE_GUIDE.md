# Student-first visual style guide

## Purpose

Every visible section should answer one physical question with one clear
interaction and one clear figure. The core reader is designed for a student on
a laptop and for discussion on a classroom projector; the same exported app
must remain usable on a phone.

## Core-figure complexity

- A core figure contains at most two panels.
- A three-panel figure is permitted only in Explore further and only when all
  three panels are necessary.
- Four-panel figures are prohibited in the core reader.
- Avoid twin y-axes in the core reader.
- Do not combine a heat map, multiple profiles, and an electrical response in
  one figure.
- One figure answers one question. If its explanation needs more than two or
  three sentences, consider splitting it.

## Typography

Use a restrained hierarchy rather than solving crowding by increasing every
font:

- base plot text: about 12--13 pt;
- axis labels: about 12--13 pt;
- titles: about 13--15 pt;
- legend text: about 10--11 pt.

Use sentence case for titles and labels. Prefer labels such as
Temperature, T (K); Bulk concentration, c∞ (cm⁻³); and
Chemical diffusivity, Dδ (cm² s⁻¹). Do not expose raw code-style labels such as
log10(D_delta / cm^2 s^-1), c_i,infinity, tau_delta, or profile boundary. When
a widget cannot render mathematics cleanly, give it a short plain-language
label and place the rendered notation in adjacent markdown.

## Lines, markers, and color

- main curve: about 1.6--1.8 pt;
- secondary curve: about 1.2--1.5 pt;
- reference line: about 1.0--1.2 pt;
- avoid 2.8--3.0 pt data lines;
- use markers sparingly;
- distinguish curves with color plus line style, marker, or direct label.

Use low-saturation colors with adequate light/dark contrast. Keep the same
physical quantity in the same color within one notebook whenever practical.

## Titles, legends, and annotations

- Keep titles short, normally under 60 characters.
- Do not put equations, full boundary conditions, or live parameter tables in
  titles.
- Do not put paragraphs inside axes or explanatory text over data.
- Place interpretation in markdown immediately below the figure.
- Use annotations only to point to a specific feature in clear empty space.
- Legends must not cover important data. Prefer direct curve labels or legends
  outside the data region when practical.

## Nyquist plots

Every Nyquist plot must:

- plot \(Z'\) horizontally and \(-Z''\) vertically;
- show units;
- use identical data scale on both axes and an equal aspect ratio;
- remain geometrically faithful when resized;
- indicate frequency direction with a few clean markers or arrows;
- never use a distorted viewport that turns a semicircle into an ellipse or
  changes a 45-degree line.

## Controls and visible-effect testing

- Show no more than three controls for one core figure.
- Put secondary and expert controls in an accordion.
- Prefer physically named presets before raw logarithmic controls.
- Remove, redesign, or move a control that produces no clear and interpretable
  visible change.
- Do not let autoscaling conceal the effect of a control.

For every core control, record the chain

\[
\text{control}\rightarrow\text{calculation}\rightarrow\text{displayed consequence}.
\]

Capture and compare the default state with low/contrasting and high/contrasting
states. If a parameter only changes scale, use fixed axes, a baseline overlay,
separate normalized and physical views, a derived scalar beside the figure, or
remove the misleading control.

## Random simulations

Do not derive a random seed from physical parameters. Use the same realization
when physical controls change so differences are causal and comparable. If a
new realization is educationally useful, provide a separate explicit action.

## Responsive browser QA

Inspect exported notebook apps at approximately:

- mobile: 390 × 844;
- laptop: 1280 × 800;
- projector: 1920 × 1080.

For every core interactive figure capture default, low/contrasting, and
high/contrasting states. Check text, controls, figures, legends, annotations,
and scrolling in the actual browser. A source review, successful export, or
matplotlib contact sheet does not count as browser visual QA.

A figure fails release if text overlaps data or other text, a legend hides an
important feature, an annotation is clipped, a control appears inert, Nyquist
geometry is distorted, or the physical point is unclear without implementation
notes.
