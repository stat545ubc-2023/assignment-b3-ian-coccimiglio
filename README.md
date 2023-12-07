
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Assignment-B3: Biological Oscillation.

## Usage

Run the Shiny App locally in the subfolder found in
/Alon_Oscillator/app.R. All of the models are compiled there.

For the original version Assignment-B3 version of this app, check out
<https://ian-fc.shinyapps.io/Alon_Oscillator/>

The updated online version for Assignment-B4 is located at
<https://ian-fc.shinyapps.io/Alon_Oscillator_Updated/>.

## New Features in Assignment B4

- Random noise checkbox for the damped oscillator, which can create
  indefinite oscillations.
- Generated second version of the repressilator with controllable
  cooperativity
- Included reactive equations for damped ODE model assessment

## Additional fixes/adaptations in Assignment B4

- Legends added to multi-element plots
- References to external work appended bottom of each model
- Exposition of input parameter functions added
- ODE Equations included to determine the parameter meanings
- Updated design and theming of figures.

## Design

I implement a series of biological oscillators as described in *Uri
Alon - An Introduction to Systems Biology* - Chapter 6. The goal is to
show how we can model these dynamical processes in R, and allow users to
moderate the parameters through an interactive GUI. Additionally, it
gives context to which conditions lead to oscillations, and allow users
to select parameters that show different dynamical behaviour.

The system is solved using the deSolve and cOde packages, which allow
for fast solving and compilation of systems of ordinary differential
equations.

### First, we build a damped oscillator.

This is a two-component system. It comes with a dynamic user interface
to edit the simulation time and the individual rate parameters. This was
modeled after Chapter 6.1 in Uri Alon’s book. In the new version of the
app, we include a “noise” checkbox that generates a vector of random,
normal noise (mean=0, sd=1) to the simulation. This causes a dynamic
shift from damped oscillations to undamped, perpetual oscillations -
similar to the dynamics in Chapter 6.2 of Alon’s book.

<div class="figure">

<img src="Alon_Oscillator/Damped.png" alt="Damped Oscillations" width="100%" />
<p class="caption">
Damped Oscillations
</p>

</div>

### Second, we showcase the Brusselator

This model has two intermediary components, X and Y, under
investigation. Despite this, the mechanism can produce stable
oscillations. This was modeled after work done by Thomas Petzoldt, here:
<https://tpetzoldt.github.io/deSolve-shiny/deSolve-shiny.html>

<div class="figure">

<img src="Alon_Oscillator/Brusselation.png" alt="Brusselator visualization" width="100%" />
<p class="caption">
Brusselator visualization
</p>

</div>

### Finally, we build two undamped 3-component oscillators

Also known as a *repressilator*. Here, we showcase a system in which
each component has negative feedback on the following component. This
has the capacity to oscillate indefinitely. We create two of these
models, one which locks the cooperativity to a set value (3), and the
second which allows free manipulation of *n*. We show that increasing
‘n’ and ‘beta’ parameters cause greater oscillations, while increasing
the alpha decay parameter leads to monotonic behaviour. The ideas for
the specific ODE system were modeled after work done by Michael Elowitz,
here: <https://biocircuits.github.io/chapters/09_repressilator.html>

<div class="figure">

<img src="Alon_Oscillator/Repression.png" alt="3 Component Oscillation" width="100%" />
<p class="caption">
3 Component Oscillation
</p>

</div>

I hope you find this work interesting, and allows you to visualize the
mathematics underlying biological oscillations!
