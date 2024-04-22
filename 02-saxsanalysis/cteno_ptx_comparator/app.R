#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(tidyverse)
library(shiny)
library(shinyWidgets)
library(shinythemes)
library(ggplot2)
library(ggpubr)

# cf https://stackoverflow.com/questions/28436842/shiny-is-it-possible-to-make-a-vertical-slider

## NOTE: this file >150 MB!
#if(exists("smoothed_data") == FALSE){
#    load("www/smoothed_JWLdata.RData")
#}

# WonB theme
theme_pubblack <- function(...){
    theme_pubr(...) +
        theme(
            # axis options
            axis.line  = element_line(color = "white"),
            axis.ticks = element_line(color = "white"),
            axis.text  = element_text(color = "white", size = 20),
            axis.title = element_text(color = "white", size = 20),
            # legend
            legend.position = "none",
            # panel
            panel.background = element_rect(color = NA,  fill = "black"),
            #panel.border = element_blank()
            # facetting
            strip.background = element_rect(fill = NA, color = "white"),
            strip.text       = element_text(color = "white"),
            # plot options
            plot.background = element_rect(color = NA,  fill = "black"),
            plot.title = element_text(color = "white")
        )
}

# COLORS
chroma = c(
    "JWL158" = "#AA7033",
    "JWL161" = "#7D91B0",
    "JWL164" = "#A8685A"
)
img_h = 100

# color font
color_font <- function(text, color){
    HTML(
        paste(
            c(
                '<font color="',
                color,
                '">',
                text,
                '</font>'
            ),
            collapse=''
        )
    )
}

# generates values for the slider
detents <- function(values){
    values = unique(values)
    c(
        "min" = min(values),
        "max" = max(values),
        "step" = (max(values) - min(values)) / (length(values)-1) %>% round(2)
    )
}

# init
#load(paste(c("www/smoothed_", "up", ".RData"), collapse=''))

# Define UI for application that draws a histogram
ui <- fluidPage(
    theme = shinytheme("cyborg"),
    setBackgroundColor("black"),
    
    titlePanel("SAXS P/T Explorer: Ctenophore Polar Lipids"),
    tags$br(),
    
    fluidRow(
        column(1,tags$img(src="Leuco_MBA.png", height = img_h)),
        column(1, offset=2, tags$h4(color_font("15°C", chroma[["JWL158"]])), align = "right"),
        column(4, 
               noUiSliderInput(
                   inputId = "temp158",
                   min = 4, 
                   max = 30,
                   step = (30-4)/24,
                   value = 15,
                   width = "250px",
                   color = chroma[["JWL158"]]
               ),
        ),
        column(1, checkboxInput("JWL158", h4(color_font("<em>Leucothea</em>", chroma[["JWL158"]])), value = TRUE))
    ),
    fluidRow(
        column(1, offset=1, tags$img(src="Bcucumis_JRW.png", height = img_h)),
        column(1, offset=1, tags$h4(color_font("0°C", chroma[["JWL161"]]), align = "right")),
        column(4, 
               noUiSliderInput(
                   inputId = "temp161",
                   min = 4, 
                   max = 60,
                   step = (60-4)/24,
                   value = 4,
                   width = "500px",
                   color = chroma[["JWL161"]]
               )
        ),
        column(1, checkboxInput("JWL161", h4(color_font("<em>Beroe</em>", chroma[["JWL161"]])), value = TRUE))
    ),
    fluidRow(
        column(1, offset=2, tags$img(src="Tjalfie_JRW.png", height = img_h)),
        column(1, tags$h4(color_font("0°C", chroma[["JWL164"]]), align = "right")),
        column(4, 
               noUiSliderInput(
                   inputId = "temp164",
                   min = 4, 
                   max = 60,
                   step = (60-4)/24,
                   value = 4,
                   width = "500px",
                   color = chroma[["JWL164"]]
               )
        ),
        column(1, checkboxInput("JWL164", h4(color_font("<em>Tjalfiella</em>", chroma[["JWL164"]])), value = TRUE))
        ),
    
    tags$br(),
    
    column(1, 
           tags$h4(color_font("1 bar", chroma[["JWL158"]]), align = "center"),
           noUiSliderInput(
               inputId = "press158", 
               min = 0, 
               max = 500,
               step = 500/24,
               value = 1,
               orientation = "vertical",
               width = "200px", 
               height = "75px",
               color = chroma[["JWL158"]]
           ),
           # dataset selector
           radioButtons(
               inputId = "dset",
               h4("Pressure Sweep"),
               list(
                   "Up"   = "up",
                   "Down" = "dn",
                   "Mean" = "av"
               ),
               selected = "up"
           ),
           textOutput(outputId = "datastatus")
    ),
    column(1, 
           tags$h4(color_font("1 bar", chroma[["JWL161"]]), align = "center"),
           noUiSliderInput(
               inputId = "press161", 
               min = 0, 
               max = 2000,
               step = 2000/24,
               value = 1,
               orientation = "vertical",
               width = "200px", 
               height = "300px",
               color = chroma[["JWL161"]]
           ),
    ),
    column(1, 
           tags$h4(color_font("400 bar", chroma[["JWL164"]]), align = "center"),
           noUiSliderInput(
               inputId = "press164", 
               min = 0, 
               max = 2000,
               step = 2000/24,
               value = 400,
               orientation = "vertical",
               width = "200px", 
               height = "300px",
               color = chroma[["JWL164"]]
           ),
    ),
    column(6,
           plotOutput("saxsProfile")
    )
)

server <- function(input, output, session) {
    
    # load selected dataset from binary
    smoothed_data = reactive({
        # ready indicator
        output$datastatus = renderText({"loading"})
        output$datastatus = renderText({"READY"})
        return(readRDS(paste(c("www/smoothed_", input$dset, ".rds"), collapse='')))
    })
    
    output$saxsProfile = renderPlot({
        smoothed_158 = smoothed_data() %>% 
            filter(input$JWL158 & samp == "JWL158") %>% 
            mutate(
                deltat = abs(temp - input$temp158),
                deltap = abs(press - input$press158)
            ) %>% 
            filter(deltat == min(deltat) & deltap == min(deltap))
        smoothed_161 = smoothed_data() %>% 
            filter(input$JWL161 & samp == "JWL161") %>% 
            mutate(
                deltat = abs(temp - input$temp161),
                deltap = abs(press - input$press161)
            ) %>% 
            filter(deltat == min(deltat) & deltap == min(deltap))
        smoothed_164 = smoothed_data() %>% 
            filter(input$JWL164 & samp == "JWL164") %>% 
            mutate(
                deltat = abs(temp - input$temp164),
                deltap = abs(press - input$press164)
            ) %>% 
            filter(deltat == min(deltat) & deltap == min(deltap))
        return(
            bind_rows(smoothed_158, smoothed_161, smoothed_164) %>% 
            #smoothed_data() %>% 
            #       filter(
            #        input$JWL158 & samp == "JWL158" &
            #          round(temp, 2) == round(input$temp158, 2) & 
            #          round(press, 2) == round(input$press158, 2) |
            #        input$JWL161 & samp == "JWL161" &
            #          round(temp, 2) == round(input$temp161, 2) & 
            #          round(press, 2) == round(input$press161, 2) |
            #        input$JWL164 & samp == "JWL164" &
            #          round(temp, 2) == round(input$temp164, 2) & 
            #          round(press, 2) == round(input$press164, 2)
            #       ) %>% 
                   ggplot(aes(x=q, y=iq, color=samp, group=samp)) +
                   geom_line(size=1.5) +
                   theme_pubblack() +
                   labs(
                       x = "q (1/Å)",
                       y = "I(q) (au)"
                   ) +
                   scale_y_log10(limits = c(min(smoothed_data()$iq), max(smoothed_data()$iq))) +
                   scale_color_manual(values = chroma)
            )
    })
    
}

shinyApp(ui, server)