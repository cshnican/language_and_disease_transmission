# Linguistic diversity as a barrier to disease spread

Sihan Chen, Antonio Benítez-Burraco

This Github repository contains the data, the code, the generated figures for the manuscript titled ''Linguistic diversity as a barrier to disease spread''.

## data
This folder contains all the dataset used in our analyses, ordered by their dataset number in Table S3 in Supplementary Information.

### **List of Datasets Used in the Analysis**

#### **Global Datasets**
- **DS1**: COVID-19 daily case data by country (2020–2023). Source: [Johns Hopkins University](https://github.com/CSSEGISandData/COVID-19/blob/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv). File: ```time_series_covid19_confirmed_global.csv```.
- **DS2**: Mean annual temperature and precipitation data by country. Source:  [Global Data Lab](https://globaldatalab.org/geos/table/). File: ```GDL-Yearly-Average-Surface-Temperature-(ºC)-data.csv``` for temperature, and ```GDL-Total-Yearly-Precipitation-(m)-data.csv``` for precipitation.
- **DS3**: Elevation by country. Source:  [Wikipedia](https://en.wikipedia.org/wiki/List_of_countries_by_average_elevation). File: ```elevations.csv```.
- **DS4**: Subnational HDI v8.1. Source:  [Global Data Lab](https://globaldatalab.org/shdi/table/). File: ```GDL-Subnational-HDI-data.csv```.
- **DS5**: Road density by country (km per surface area). Source:  [World Bank](https://prosperitydata360.worldbank.org/en/indicator/WEF+TTDI+ROADDENS). File: ```WEF-TTDI.csv```.
- **DS6**: Population density by country. Source:  [World Bank](https://data.worldbank.org/indicator/EN.POP.DNST). File: ```API_EN.POP.DNST_DS2_en_csv_v2_89.csv```
- **DS7**: Population by country (2022). Source:  [World Bank](https://data.worldbank.org/indicator/SP.POP.TOTL). File: ```populations.csv```.
- Language diversity index (LDI): Source: [Ethnologue](https://www.ethnologue.com/statistics/). File: ```ethonologue_diversity.csv```.

#### **U.S. Datasets**
- **DS8**: "Languages spoken at home" (ACS 2023 5-year estimates). Source:  [U.S. Census Bureau](https://data.census.gov/table?q=cp02&g=010XX00US$0500000). File: ```ACSCP5Y2023.CP02_2025-02-06T154350/ACSCP5Y2023.CP02-Data.csv```. See entries labelled as: ```CP02_2023_113E, CP02_2023_114E, CP02_2023_115E, CP02_2023_116E, CP02_2023_117E, CP02_2023_118E, CP02_2023_119E, CP02_2023_120E,, CP02_2023_121E, CP02_2023_122E, CP02_2023_123E```.
- **DS9**: Mean annual temperature by U.S. county (2022). Source:  [NOAA](https://www.ncei.noaa.gov/access/monitoring/climate-at-a-glance/county/mapping/110/tavg/202212/12/value). File: ```county_temp.csv``` for the Lower 48, ```county_temp_AK.csv``` for Alaska, ```county_temp_HI.csv``` for Hawaii.
- **DS10**: Mean annual precipitation by U.S. county (2022). Source:  [NOAA](https://www.ncei.noaa.gov/access/monitoring/climate-at-a-glance/county/mapping/110/pcp/2022/12/value). File: ```county_precip.csv``` for the Lower 48, ```county_precip_AK.csv``` for Alaska, ```county_precip_HI.csv``` for Hawaii.
- **DS11**: Latitude and area of each U.S. county. Source:  [U.S. Census Bureau](https://www2.census.gov/geo/docs/maps-data/data/gazetteer/2024_Gazetteer/2024_Gaz_counties_national.zip). File: ```2024_Gaz_counties_national.txt```.
- **DS12**: "Education attainment" data (ACS 2023 5-year estimates). Source:  [U.S. Census Bureau](https://data.census.gov/table/ACSCP5Y2023.CP02?q=CP02:+Comparative+Social+Characteristics+in+the+United+States). File: ```ACSCP5Y2023.CP02_2025-02-06T154350/ACSCP5Y2023.CP02-Data.csv```. See entries labeled as: ```CP02_2023_060E, CP02_2023_061E, CP02_2023_062E, CP02_2023_063E, CP02_2023_064E, CP02_2023_065E, CP02_2023_066E```.
- **DS13**: "Income in the past 12 months" (ACS 2023 5-year estimates). Source:  [U.S. Census Bureau](https://data.census.gov/table?q=S1901:+Income+in+the+Past+12+Months+(in+2023+Inflation-Adjusted+Dollars)&y=2023&d=ACS+5-Year+Estimates+Subject+Tables). File: ```ACSST5Y2023.S1901_2025-02-06T234736/ACSST5Y2023.S1901-Data.csv```.
- **DS14**: Life expectancy by county (2010-2015). Source:  [CDC](https://www.cdc.gov/nchs/data-visualization/life-expectancy/index.html). File: ```U.S._Life_Expectancy_at_Birth_by_State_and_Census_Tract_-_2010-2015.csv```.
- **DS15**: Population by county (2023). Source:  [U.S. Census Bureau](https://www2.census.gov/programs-surveys/popest/datasets/2020-2023/counties/totals/co-est2023-alldata.csv). File: ```co-est2023-alldata.csv```.
- **DS16**: Road length by county, compiled by Erin Davis. Source:  [U.S. Census Bureau](https://gist.github.com/erdavis1/c275afbcfffbb6bf01ade657ad0ddd65). File: ```road_length_by_county.csv```.
- **DS17**: List of the busiest passenger flight routes. Source:  [Wikipedia](https://en.wikipedia.org/wiki/List_of_busiest_passenger_flight_routes). 

## scripts

- ```Global_analysis.R```: script for the global analysis.
- ```US_analysis.R```: script for the US analysis.
- ```run_analysis.m```: script for the simulation.
- ```plot_halfpoint.R```: script for plotting halfpoints in Fig. 3.
- ```trudgill_network_spatiotemporal_ODE_solver.m```: scripts for solving the ODE
- ```trudgill_network_spatiotemporal_ODE.m```: scripts laying out the ODE


## figs

- ```beta_us.pdf```: Figure S2.
- ```PDE_illu.pdf```: Figure 2.
- ```reg_us.pdf```: Figure 1b.
- ```reg.pdf```: Figure 1a.
- ```scenarios_1-6_hp.pdf```: Figure 3.
- ```scenarios_12.pdf```: Figure S3.
- ```scenarios_34.pdf```: Figure S4.
- ```scenarios_56.pdf```: Figure S5.
- ```transmission_rate```: Figure S1.

## output_table

- ```half_points_scenario1-6.txt```: Half points for Scenarios 1-6, data used for Fig. 3.
- ```scenarios1.txt``` - ```scenarios6.txt```: Complete disease transmission dynamics in each scenario.

