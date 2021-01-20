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


## Server components are defined here

## The global definitions are in global.R!

## The data is in the following variables:
## data$r real matrix
## data$f list of factors
## data$s subsets extracted from the first factor
## The state that may vary is collected to current as follows:
##  current$s the current subset
##  current$w information about the current projection
##  current$v PCA eigenvalues
## In your application the current may containt something totally
## different.

summarizetiles <- function(current) {
    R <- current$R
    C <- current$C
    tiles <- current$tiles
    a <- data.frame(i=0,name="current",sizeR=length(R),sizeC=length(C),
                    recall=1,precision=1,jaccard=1,H1=FALSE,H2=FALSE)

    if(length(tiles)>0) {
        a <- rbind(a,
                   data.frame(i=1:length(tiles),
                              name=names(tiles),
                              sizeR=sapply(tiles,function(x) length(x$R)),
                              sizeC=sapply(tiles,function(x) length(x$C)),
                              recall=round(sapply(tiles,function(x) length(intersect(x$R,R))/length(x$R)),digits=3),
                              precision=round(sapply(tiles,function(x) length(intersect(x$R,R))/length(R)),digits=3),
                              jaccard=round(sapply(tiles,function(x) length(intersect(x$R,R))/length(union(x$R,R))),digits=3),
                              H1=sapply(current$tiles,function(x) x$H1),
                              H2=sapply(current$tiles,function(x) x$H2)))
    }
    a
}


shinyServer(function(input,output,session) {
    ## This function updates plots etc. This should be fast, i.e.,
    ## no computations.
    ## We call it immediately after observing an input event which
    ## should update plots.
    updateall <- function() {
      withProgress(message="updating view...",value=0,{
        output$plot <<- renderPlot({plotdata(data,current)})
        output$plotattr <<- renderPlot({plotattr(data,current,
                                                 factors=current$factors)})
        output$pcplot <<- renderPlot({pcplot(data,current,current$pc)})
        output$seldescr <<- renderText(sprintf("%d of %d rows selected",
                                       length(current$R),
                                       dim(data$r)[1]))
        output$seldescr2 <<- renderText(sprintf("%d of %d columns selected",
                                       length(current$C),
                                       dim(data$data)[2]))
        output$datadescr <<- renderText(data$descr)
        output$tiledescr <<- renderText(
            sprintf("H1 %d tiles, H2 %d tiles",
                    length(ls(current$H1$status()$Tl)),
                    length(ls(current$H2$status()$Tl))))
        ## output$tiles <<- renderDT(summarizetiles(current))
        output$tiles <<- DT::renderDataTable(
            summarizetiles(current),
            selection=list(mode="single",target="cell"),
            rownames=FALSE
        )
        updateCheckboxGroupInput(session,"checkGroup",choices = makecolumnchoises(),selected=current$C)
        ## print(current$C)
      })
    }

    updatetiling <- function() {

        fH1 <- sapply(current$tiles,function(x) x$H1)
        fH2 <- sapply(current$tiles,function(x) x$H2)
        H1 <- tiling(dim(data$data)[1],dim(data$data)[2])
        for(i in which(fH1 & fH2)) {
          H1$addtile(R=current$tiles[[i]]$R,C=current$tiles[[i]]$C)
        }
        H2 <- H1$copy()
        for(i in which(fH1 & !fH2)) {
          H1$addtile(R=current$tiles[[i]]$R,C=current$tiles[[i]]$C)
        }
        for(i in which(!fH1 & fH2)) {
          H2$addtile(R=current$tiles[[i]]$R,C=current$tiles[[i]]$C)
        }

      current$H1 <<- H1
      current$H2 <<- H2
    }

    observeEvent(input$tiles_cells_selected,{
        ij <- input$tiles_cells_selected
        if(dim(ij)[1]>0) {
            i <- ij[1,1]-1
            j <- ij[1,2]
            if(i>0 && j %in% c(2,3,7,8)) {
                if(j==2) { # R
                    current$R <<- union(current$R,current$tiles[[i]]$R)
                } else if(j==3) { #C
                    current$C <<- union(current$C,current$tiles[[i]]$C)
                } else if(j==7) { # H1
                    current$tiles[[i]]$H1 <<- !(current$tiles[[i]]$H1)
                    updatetiling()
                } else if(j==8) { # H2
                    current$tiles[[i]]$H2 <<- !(current$tiles[[i]]$H2)
                    updatetiling()
                }


                updateall()
            }
        }
    })

    observeEvent(input$sel1,{
        s <- input$sel1
        xy <- data$r %*% current$w[,1:2]
        a <- which(s$xmin<=xy[,1] & xy[,1]<=s$xmax &
                       s$ymin<=xy[,2] & xy[,2]<=s$ymax)
        current$R <<- if(input$rowmode=="add") union(current$R,a) else setdiff(current$R,a)
        updateall()
    })

    observeEvent(input$selattr,{
        s <- input$selattr
        xy <- data$raux %*% current$w[,1:2]
        a <- which(s$xmin<=xy[,1] & xy[,1]<=s$xmax &
                       s$ymin<=xy[,2] & xy[,2]<=s$ymax)
        current$C <<- if(input$colmode=="add") union(current$C,a) else setdiff(current$C,a)
        updateall()
    })

    observeEvent(input$sel2,{
        s <- input$sel2
        i <- min(length(current$pc),max(1,ceiling(s$ymin)))
        j <- min(length(current$pc),max(1,floor(s$ymax)))
        if(i<=j) {
            a <- current$pc[i:j]
            current$C <<- if(input$pcsmode=="add") union(current$C,a) else setdiff(current$C,a)
        }
        updateall()
    })

    observeEvent(input$clearselection,{
        current$R <<- c()
        updateall()
    })

    observeEvent(input$clearselection2,{
        current$C <<- c()
        updateall()
    })


    observeEvent(input$reverseselection,{
        current$R <<- setdiff(1:dim(data$r)[1],current$R)
        updateall()
    })

    observeEvent(input$reverseselection2,{
        current$C <<- setdiff(1:dim(data$data)[2],current$C)
        updateall()
    })

    observeEvent(input$togglefactors,{
        current$factors <<- !current$factors
        updateall()
    })

    observeEvent(input$recompute,{
        cat("server.R: computing PCA.\n")
        withProgress(message="computing...",value=0,{
          a <- dopca(if(length(current$R)>1) data$r[current$R,] else data$r)
        })
        current$label <<- "PC"
        current$w <<- a$w
        current$v <<- a$v
        updateall()
    })

    observeEvent(input$recomputet,{
        cat("server.R: computing tiling projection.\n")
        ii <- 1:dim(data$r)[2]
        withProgress(message="computing...",value=0,{
          a <- findproj2(current$H1$cov(data$data)[ii,ii],
                         current$H2$cov(data$data)[ii,ii])
        })
        current$label <<- "T"
        current$w <<- normw(a$w)
        current$v <<- a$v
        updateall()
    })

    observeEvent(input$addtile,{
        current$tiles[[sprintf("user%02d",
                               current$count)]] <<- list(R=current$R,
                                                         C=current$C,
                                                         H1=TRUE,
                                                         H2=TRUE)
	updatetiling()
        current$count <<- current$count+1
        updateall()
    })

    observeEvent(input$datafile,{
        data <<- getdata(input$datafile)
        current <<- initcurrent(data)
        updateSliderInput(session,"pccols",
                          min=1,max=dim(data$r)[2],
                          value=min(10,dim(data$r)[2]))
        updateall()
    })

    observeEvent(input$updatetilings,{
        updatetilings()
        updateall()
    })

    observeEvent(input$pccols,{
        current$pc <<- pcorder(data,current,current$pcmode,input$pccols)
        updateall()
    })

    observeEvent(input$pcmodeX1,{
        current$pc <<- pcorder(data,current,"X1",input$pccols)
        current$pcmode <<- "X1"
        updateall()
    })

    observeEvent(input$pcmodesel,{
        current$pc <<- pcorder(data,current,"selection",input$pccols)
        current$pcmode <<- "selection"
        updateall()
    })

    observeEvent(input$datadump,{
        saveRDS(data,file="DATA.rds")
        cat("Dumped DATA.rds.\n")
    })

    observeEvent(input$currentdump,{
        saveRDS(current,file="CURRENT.rds")
        cat("Dumped CURRENT.rds.\n")
    })

    observeEvent(input$update,{
      current$C <<- as.numeric(input$checkGroup)
      updateall()
    })
})
