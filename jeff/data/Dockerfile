FROM lc_app_upper_base_04jul2020
WORKDIR /home/shiny
COPY jeffApp.R .
COPY tbl.39analytes.1157x6.RData .
COPY installPrerequisites.R .
RUN R -f installPrerequisites.R
COPY pkgs/igvShiny_0.99.91.tar.gz .
COPY pkgs/shinyModules_1.0.13.tar.gz .
RUN R CMD INSTALL igvShiny_0.99.91.tar.gz
RUN R CMD INSTALL shinyModules_1.0.13.tar.gz



