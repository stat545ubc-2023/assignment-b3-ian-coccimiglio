#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(cOde)
library(deSolve)

# Feature 1, Creating a damped biological oscillator. This compiles an object which
# simulates a biological system that tends towards damped oscillation.
# Users can change the INDIVIDUAL parameters and the simulation time.
# Outcomes are changes in the phase-space diagram and the time-concentration plot.
osc2 <- cOde::funC(
    c(
        X = "-b1*Y - a1*X",
        Y = "b2*X - a2*Y"
    ), modelname = "Alon_damped", compile=TRUE
)

## Feature 2, Creating a 3-part repressilator. This compiles an object which
# simulates a 3 component system that tends towards infinite oscillations.
# Users are able to change the global degradation, creation, and concentration constants.
# Additionally, they can change the simulation length.
# As above, this produces a phase-space diagram and the time-concentration plot.
osc3 <- cOde::funC(
    c(
        # 1 and 3 are rate constants for the purposes of this demo.
        X = "beta/(1+(Z/k)^3) - g2*X",
        Y = "beta/(1+(X/k)^3) - g2*Y",
        Z = "beta/(1+(Y/k)^3) - g2*Z"
    ),
    modelname = "Osc_3", compile=TRUE
)

## Feature 3, Creating the brusselator This compiles an object which
# simulates a 2 component system that tends towards infinite oscillations.
# Users are able to manipulate the different inflows and outflow rates,
# as well as concentrations of products A, and B
# As above, this produces a phase-space diagram and the time-concentration plot.
brussC <- cOde::funC(
    c(
        X = "k1*A - k2*B*X + k3*(X^2)*Y - k4*X",
        Y = "k2*B*X - k3*(X^2)*Y"
    ), modelname = "brusselC", compile=TRUE
)


# Define UI for application that draws a histogram
ui <- fluidPage(
    titlePanel("Biological Oscillations"),
    tabsetPanel(
        tabPanel("Damped Biological Oscillator",
            sidebarLayout(
                sidebarPanel(
                    sliderInput("alpha1", "Degradation X (\u03B11)", min = 1, max = 20, value = 1),
                    sliderInput("beta1", "Formation X (\u03B21)", min = 1, max = 20, value = 8),
                    sliderInput("alpha2", "Degradation Y (\u03B12)", min = 1, max = 20, value = 1),
                    sliderInput("beta2", "Formation Y (\u03B22)", min = 1, max = 20, value = 8),
                    sliderInput("len", "Simulation Time", min = 1, max = 50, value = 6)
                ),

                mainPanel(
                         plotOutput("dampOsc"),
                         imageOutput("damp"),
                )
            )
        ),

        tabPanel("Repressilator",
                 sidebarLayout(
                     sidebarPanel(
                         sliderInput("Beta", "Protein Generation Rate (\u03B2)", min = 0, max = 20, value = 1),
                         sliderInput("k", "Protein Concentration Factor", min = 1, max = 20, value = 1),
                         sliderInput("Gamma", "Protein Degradation Rate (\u03B3)", min = 0, max = 5, value = 0.1, step=0.1),
                         sliderInput("len_2", "Simulation Length", min = 1, max = 1000, value = 300)
                     ),

                     # Show a plot of the generated distribution
                     mainPanel(
                         plotOutput("Osc3"),
                         imageOutput("repress")
                     )
                 )
        ),


        tabPanel("Brusselator",
             sidebarLayout(
                 sidebarPanel(
                     sliderInput("A", "Autocatalytic A", min = 1, max = 20, value = 1),
                     sliderInput("B", "Autocatalytic B", min = 1, max = 20, value = 3),
                     sliderInput("k1", "Catalytic Rate A (k1)", min = 1, max = 20, value = 1),
                     sliderInput("k2", "Catalytic Rate B (k2)", min = 1, max = 20, value = 1),
                     sliderInput("k3", "XY binding rate (k3)", min = 1, max = 20, value = 1),
                     sliderInput("k4", "X degradation rate (k4)", min = 1, max = 20, value = 1),
                     sliderInput("len_3", "Simulation Time", min = 1, max = 200, value = 50)
                 ),

                 # Show a plot of the generated distribution
                 mainPanel(
                     plotOutput("Brussel"),
                     imageOutput("brussel")
                 )
             )
        )

    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$damp <- renderImage({
        list(src = "Damped_Oscillator.png",
        width = 200,
        height = 200)
    }, deleteFile = FALSE)

    output$repress <- renderImage({
        list(src = "Repressilator.png",
             width = 200,
             height = 250)
    }, deleteFile = FALSE)

    output$brussel <- renderImage({
        list(src = "Brusselator.png",
             width = 250,
             height = 250)
    }, deleteFile = FALSE)

    output$dampOsc <- renderPlot({
        par(mfrow=c(1,2))
        y0 <- c(X=1,Y=1)
        k_list <- c(b1=input$beta1, b2=input$beta2, a1=input$alpha1, a2=input$alpha2)
        step=0.01
        damp_osc <- odeC(y=y0, times=seq(0,input$len,step), func=osc2, parms = k_list)
        plot(damp_osc[,1], damp_osc[,2], type='l',
             xlab = "Time",
             ylab = "Concentration")
        lines(damp_osc[,1], damp_osc[,3], col='red', type='l')
        plot(damp_osc[,2], damp_osc[,3], type='l',
             xlab = "Concentration X",
             ylab = "Concentration Y"
             )
    }
    )

    output$Osc3 <- renderPlot({
        par(mfrow=c(1,2))
        Z_Conc = 3
        Y_Conc = 2
        X_Conc = 1
        y0_2 <- c(X=X_Conc,Y=Y_Conc,Z=Z_Conc)
        k_list_2 <- c(beta=input$Beta, k=input$k, g2=input$Gamma)
        step=0.2
        osc_3 <- odeC(y=y0_2, times=seq(0,input$len_2,step), func=osc3, parms = k_list_2)
        plot(osc_3[,1], osc_3[,2], col='black', type='l',
             xlab = "Time",
             ylab = "Concentration")
        lines(osc_3[,1], osc_3[,3], col='red')
        lines(osc_3[,1], osc_3[,4], col='blue')
        plot(osc_3[,2], osc_3[,3], col='black', type='l',
             xlab = "Concentration X",
             ylab = "Concentration Y"
             )
    }
    )

    output$Brussel <- renderPlot({
        par(mfrow=c(2,1))
        X_Conc_Bruss <- 2
        Y_Conc_Bruss <- 2
        y0_3 <- c(X=X_Conc_Bruss,Y=Y_Conc_Bruss)
        k_list3 <- c(A=input$A, B=input$B, k1=input$k1, k2=input$k2, k3=input$k3, k4=input$k4)
        step=0.1
        ode_bruss <- odeC(y=y0_3, times=seq(0,input$len_3,step), func=brussC, parms = k_list3)
        plot(ode_bruss[,1], ode_bruss[,2],col='black', type='l',ylim=c(0,8),
             xlab = "Time",
             ylab = "Concentration")
        lines(ode_bruss[,1], ode_bruss[,3],col='blue')
        plot(ode_bruss[,2], ode_bruss[,3],type='l',
             xlab = "Concentration X",
             ylab = "Concentration Y")
    }
    )
}

# Run the application
shinyApp(ui = ui, server = server)
