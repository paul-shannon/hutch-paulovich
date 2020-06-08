FROM rocker/shiny
WORKDIR /home/shiny
COPY tbl.numa.RData .
COPY firstApp.R .
COPY installPrerequisites.R .
RUN R -f ./installPrerequisites.R .
CMD ["R", "-e", "shiny::runApp('/home/shiny/lcGenes.R', port=3838, host='0.0.0.0')"]





