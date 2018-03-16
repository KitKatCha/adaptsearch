#load packages
library(shiny)
library(shinyWidgets)
library(xcms)
library(plotly)
library(stringr)


load("/srv/shiny-server/data/inputdata.dat")
raw_files <- basename(rownames(xdata@phenoData@data))

# identifier_type name because of filename and not hid
gx_get(raw_files, identifier_type='name')

for (file in raw_files){
  system(sprintf("ln -s /import/%s %s", file, file))
}

xdata <- adjustRtime(xdata, param = ObiwarpParam(binSize = 0.6))
adjusted <- hasAdjustedRtime(xdata)


ui <- bootstrapPage(
	conditionalPanel(
		condition = "output.hide_panel",
		switchInput(
			inputId = "adjustedTime",
			label = strong("Adjusted Time"),
			value = TRUE
		)
	),
	fluidRow(
		column(
			6,
			div(style = "border-style: solid; border-width: thin; border-color: #000000",
				plotlyOutput('TIC')
			)
		),
		column(
			6,
			div(style = "border-style: solid; border-width: thin;border-color: #000000",
				plotlyOutput('TIC_chromato')
			)
		)
	),
	fluidRow(
		column(
			6,
			div(style = "border-style: solid; border-width: thin; border-color: #000000",
				plotlyOutput('BPC')
			)
		),
		column(
			6,
			div(style = "border-style: solid; border-width: thin; border-color: #000000",
				plotlyOutput('BPC_chromato')
			)
		)
	)
)

server <- function(input, output){

	output$hide_panel <- eventReactive( adjusted, TRUE, ignoreInit = TRUE)
	outputOptions(output, "hide_panel", suspendWhenHidden = FALSE)
	
	output$TIC <- renderPlotly({

                # According to Adjusted Time
		if (input$adjustedTime) {
			title <- "TIC adjusted"
			rtime <- split(rtime(xdata, adjusted = TRUE)/60, f = fromFile(xdata))
                } else {
			title <- "TIC with tic()"
			rtime <- split(rtime(xdata, adjusted = FALSE)/60, f = fromFile(xdata))
		}

		#function tic of MsnBase faster than chromatogram of XCMS (for me, to check)
		intensity <- split(tic(xdata), f = fromFile(xdata))

		chromato <- plot_ly(source='alignmentChromato', type='scatter', mode='markers') %>% 
			layout(title=title, xaxis=list(title='Retention time (min)'), yaxis=list(title='Intensity'), showlegend=TRUE) %>% 
			config(scrollZoom=TRUE, showLink=TRUE, displaylogo=FALSE, 
				modeBarButtons=list(list('toImage', 'zoom2d', 'select2d', 'pan2d', 'autoScale2d', 'resetScale2d')))
		
		if(is.null(raw_files)) return(chromato)
		else if(!length(raw_files)) return(chromato)
		
		for(i in 1:length(raw_files)) {
			index <- which(basename(rownames(phenoData(xdata))) == raw_files[i])
			chromato <- chromato %>% add_lines(
				x=rtime[[index]], y=intensity[[index]], name=basename(raw_files[i]), hoverinfo='text', 
				text=paste('Intensity: ', round(intensity[[index]]), '<br />Retention Time: ', round(rtime[[index]], digits=2))
			)
		}

		return(chromato)
	})

        output$TIC_chromato <- renderPlotly({

                # According to Adjusted Time
                if (input$adjustedTime) {
                        title <- "TIC adjusted with chromatogram()"
                        rtime <- chromatogram(xdata, aggregationFun = 'sum', adjustedRtime = TRUE)
                } else {
                        title <- "TIC with chromatogram()"
                        rtime <- chromatogram(xdata, aggregationFun = 'sum', adjustedRtime = FALSE)
                }

                chromato <- plot_ly(source='alignmentChromato', type='scatter', mode='markers') %>%
                        layout(title=title, xaxis=list(title='Retention time (min)'), yaxis=list(title='Intensity'), showlegend=TRUE) %>%
                        config(scrollZoom=TRUE, showLink=TRUE, displaylogo=FALSE,
                                modeBarButtons=list(list('toImage', 'zoom2d', 'select2d', 'pan2d', 'autoScale2d', 'resetScale2d')))

                if(is.null(raw_files)) return(chromato)
                else if(!length(raw_files)) return(chromato)

                for(i in 1:length(raw_files)) {
                        index <- which(basename(rownames(phenoData(xdata))) == raw_files[i])
                        chromato <- chromato %>% add_lines(
                                x=rtime[[index]]@rtime/60, y=rtime[[index]]@intensity, name=basename(raw_files[i]), hoverinfo='text',
                                text=paste('Intensity: ', round(rtime[[index]]@intensity), '<br />Retention Time: ', round(rtime[[index]]@rtime/60, digits=2))
                        )
                }

                return(chromato)
        })


	
	output$BPC <- renderPlotly({

                # According to Adjusted Time
                if (input$adjustedTime) {
                        title <- "BPC adjusted"
                        rtime <- split(rtime(xdata, adjusted = TRUE)/60, f = fromFile(xdata))
                } else {
                        title <- "BPC with bpi()"
                        rtime <- split(rtime(xdata, adjusted = FALSE)/60, f = fromFile(xdata))
                }

                #function tic of MsnBase faster than chromatogram of XCMS (for me, to check)
                intensity <- split(bpi(xdata, initial = FALSE), f = fromFile(xdata))

		chromato <- plot_ly(source='alignmentChromato', type='scatter', mode='markers') %>% 
			layout(title=title, xaxis=list(title='Retention time (min)'), yaxis=list(title='Intensity'), showlegend=TRUE) %>% 
			config(scrollZoom=TRUE, showLink=TRUE, displaylogo=FALSE, 
				modeBarButtons=list(list('toImage', 'zoom2d', 'select2d', 'pan2d', 'autoScale2d', 'resetScale2d')))

		if(is.null(raw_files)) return(chromato)
		else if(!length(raw_files)) return(chromato)

                for(i in 1:length(raw_files)) {
                        index <- which(basename(rownames(phenoData(xdata))) == raw_files[i])
                        chromato <- chromato %>% add_lines(
                                x=rtime[[index]], y=intensity[[index]], name=basename(raw_files[i]), hoverinfo='text',
                                text=paste('Intensity: ', round(intensity[[index]]), '<br />Retention Time: ', round(rtime[[index]], digits=2))
                        )
                }

		return(chromato)
	})

        output$BPC_chromato <- renderPlotly({

                # According to Adjusted Time
                if (input$adjustedTime) {
                        title <- "BPC adjusted with chromatogram()"
	                rtime <- chromatogram(xdata, aggregationFun = 'max', adjustedRtime = TRUE)
                } else {
                        title <- "BPC with chromatogram()"
	                rtime <- chromatogram(xdata, aggregationFun = 'max', adjustedRtime = FALSE)
                }

                chromato <- plot_ly(source='alignmentChromato', type='scatter', mode='markers') %>%
                        layout(title=title, xaxis=list(title='Retention time (min)'), yaxis=list(title='Intensity'), showlegend=TRUE) %>%
                        config(scrollZoom=TRUE, showLink=TRUE, displaylogo=FALSE,
                                modeBarButtons=list(list('toImage', 'zoom2d', 'select2d', 'pan2d', 'autoScale2d', 'resetScale2d')))

                if(is.null(raw_files)) return(chromato)
                else if(!length(raw_files)) return(chromato)

                for(i in 1:length(raw_files)) {
                        index <- which(basename(rownames(phenoData(xdata))) == raw_files[i])
                        chromato <- chromato %>% add_lines(
                                x=rtime[[index]]@rtime/60, y=rtime[[index]]@intensity, name=basename(raw_files[i]), hoverinfo='text',
                                text=paste('Intensity: ', round(rtime[[index]]@intensity), '<br />Retention Time: ', round(rtime[[index]]@rtime/60, digits=2))
                        )
                }

                return(chromato)
        })


}

shinyApp(ui, server)
