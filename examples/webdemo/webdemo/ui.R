## corand
##
## Copyright 2019 Kai Puolamaki <kai.puolamaki@iki.fi>
##
## Permission is hereby granted, free of charge, to any person obtaining
## a copy of this software and associated documentation files (the
## "Software"), to deal in the Software without restriction, including
## without limitation the rights to use, copy, modify, merge, publish,
## distribute, sublicense, and/or sell copies of the Software, and to
## permit persons to whom the Software is furnished to do so, subject to
## the following conditions:
##
## The above copyright notice and this permission notice shall be
## included in all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
## EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
## MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
## NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
## LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
## OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
## WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## UI components are defined here

## You can define page layouts using most HTML widgets here.
## See https://shiny.rstudio.com/gallery/

library("shiny")


shinyUI(fluidPage(
  titlePanel("CORAND demo"),
  sidebarLayout(
      sidebarPanel(
          HTML("<strong>select rows</strong>"),
          htmlOutput("seldescr"),
          actionButton("clearselection",label="clear current selection"),
          actionButton("reverseselection",label="reverse current selection"),
          selectInput("rowmode",label="selection mode",
                      choices=c("add","remove")),
          hr(),
          HTML("<strong>select columns</strong>"),
          htmlOutput("seldescr2"),
          actionButton("clearselection2",label="clear current selection"),
          actionButton("reverseselection2",label="reverse current selection"),
          actionButton("togglefactors",label="toggle factors"),
          selectInput("colmode",label="selection mode",
                      choices=c("add","remove")),
          hr(),
          HTML("<strong>recompute projection</strong>"),
          actionButton("recompute",label="PCA projection on selection"),
          actionButton("recomputet",label="tiling projection"),
          hr(),
          HTML("<strong>parallel coordinates</strong>"),
          sliderInput("pccols","columns:",
                      min=1,max=dim(data$r)[2],
                      value=min(10,dim(data$r)[2])),
          actionButton("pcmodeX1",label="order by X1"),
          actionButton("pcmodesel",label="order by selection"),
          selectInput("pcsmode",label="selection mode",
                      choices=c("add","remove")),
          hr(),
          HTML("<strong>tiling</strong>"),
          ##actionButton("updatetilings",label="update tilings"),
          htmlOutput("tiledescr"),
          hr(),
          HTML("<strong>tiles</strong>"),
          actionButton("addtile",label="add current tile"),
          hr(),
          selectInput("datafile",label="datafile",choices=datafiles),
          htmlOutput("datadescr"),
          actionButton("datadump",label="dump DATA.rds"),
          actionButton("currentdump",label="dump CURRENT.rds"),
          hr(),
          HTML("<p>See <a href=\"https://github.com/edahelsinki/corand\">https://github.com/edahelsinki/corand</a> for more information.</p>")
    ),
    mainPanel(
        plotOutput("plot",brush=brushOpts(id="sel1"),click="plot_click"),
        plotOutput("plotattr",brush=brushOpts(id="selattr"),click="plot_click"),
        plotOutput("pcplot",brush=brushOpts(id="sel2"),click="plot_click"),
        DT::dataTableOutput("tiles"),
        checkboxGroupInput("checkGroup",
          label = h3("columns"),
          choices = makecolumnchoises(),
          selected = current$C),
        actionButton("update",label="update column selection")
    )
  )
))
