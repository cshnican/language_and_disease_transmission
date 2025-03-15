library(tidyverse)
library(ggplot2)

county_regex <- "(.*?)(?=\\s+(County|Borough|Municipality|Census Area|Parish|Planning Region|city|City and Borough)$)"

county_map_data <- usmap::us_map(regions='counties')

cosine_similarity <- function(x, y) {
  sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))
}

population <- read.csv('../data/co-est2023-alldata.csv') %>%
  filter(COUNTY != 0) %>%
  dplyr::select(STNAME, CTYNAME, POPESTIMATE2022) %>%
  mutate(Admin2 = str_extract(CTYNAME, county_regex))
  

# NaNs are mainly from out of state populations
# exclude Puerto Rico, city data from VA, Dept of Correction data in each state
# exclude certain counties in UT where no cases has been reported

timeseries_data <- read.csv('../data/time_series_covid19_confirmed_US.csv') %>% 
  pivot_longer(starts_with('X'), names_to = 'date', values_to = 'cumulative_cases') %>%
  mutate(date = str_sub(date, 2, -1),
         date = as.Date(date, '%m.%d.%y')) %>%
  separate(date, into=c('year', 'month', 'day'), remove=FALSE) %>%
  group_by(Province_State, Admin2, Combined_Key, date, month, day, year) %>%
  summarize(cumulative_cases = sum(cumulative_cases)) %>%
  left_join(population, by=c('Province_State' = 'STNAME', 'Admin2')) %>% ungroup() %>%
  group_by(Province_State, Admin2, Combined_Key) %>%
  mutate(IN = c(0, diff(cumulative_cases)),
         S = POPESTIMATE2022 - cumulative_cases,
         beta = (-log(1 - IN/(cumulative_cases + 1e-30) + 1e-30))/(1*S/POPESTIMATE2022)) %>%
  summarize(mean_beta = mean(beta, na.rm = TRUE)) %>% 
  filter(!is.nan(mean_beta),
        mean_beta > 0) 

# social demographic data: https://data.census.gov/table?q=cp02&g=010XX00US$0500000
# also remove Puerto Rico
cp02 <- read.csv('../data/ACSCP5Y2023.CP02_2025-02-06T154350/ACSCP5Y2023.CP02-Data.csv') 

language_diversity <- cp02 %>%
  filter(!grepl(', Puerto Rico', NAME)) %>%
  select(NAME, starts_with('CP02_2023_11'), CP02_2023_120E, CP02_2023_121E, CP02_2023_122E, CP02_2023_123E) %>%
  select(-CP02_2023_110E, -CP02_2023_111E, -CP02_2023_112E) %>%
  rename(
    English_only = CP02_2023_113E,
    `Language other than English` = CP02_2023_114E,
    `Language other than English!!Speak English less than "very well"` = CP02_2023_115E,
    Spanish = CP02_2023_116E,
    `Spanish!!Speak English less than "very well"` = CP02_2023_117E,
    `Other Indo-European languages` = CP02_2023_118E,
    `Other Indo-European languages!!Speak English less than "very well"` = CP02_2023_119E,
    `Asian and Pacific Islander languages` = CP02_2023_120E,
    `Asian and Pacific Islander languages!!Speak English less than "very well"` = CP02_2023_121E,
    `Other languages` = CP02_2023_122E,
    `Other languages!!Speak English less than "very well"` = CP02_2023_123E
  )


### education
education_profile <- cp02 %>% 
  filter(!grepl(', Puerto Rico', NAME)) %>%
  select(NAME, CP02_2023_060E, CP02_2023_061E, CP02_2023_062E, CP02_2023_063E, CP02_2023_064E, CP02_2023_065E, CP02_2023_066E) %>%
  mutate(across(-1, function(x) as.numeric(x)/100)) 

counties <- education_profile$NAME
education_profile_m <- education_profile %>% select(-NAME) %>% as.matrix()

education_similarity <- outer(1:nrow(education_profile), 1:nrow(education_profile), 
                       Vectorize(function(i, j) cosine_similarity(education_profile_m[i, ], education_profile_m[j, ]))) %>%
  as.table() %>% as.data.frame() %>%
  cbind(expand_grid(county2 = counties, county1 = counties)) %>%
  mutate(county_pair = pmap_chr(list(county1, county2), 
                                 ~ paste(sort(c(..1, ..2)), collapse = "+"))) %>%
  filter(county1 != county2) %>%
  select(-county1, -county2, -Var1, -Var2) %>% distinct() %>%
  separate(county_pair, into=c('county.name1', 'county.name2'), sep='\\+', remove=FALSE) %>%
  separate(county.name1, into=c('county1', 'state1'), sep=', ', remove=FALSE) %>%
  separate(county.name2, into=c('county2', 'state2'), sep=', ', remove=FALSE) 

### income
income_profile <- read.csv('../data/ACSST5Y2023.S1901_2025-02-06T234736/ACSST5Y2023.S1901-Data.csv') %>%
  filter(!grepl(', Puerto Rico', NAME),
         grepl(',', NAME)) %>%
  select(NAME, S1901_C01_002E, S1901_C01_003E, S1901_C01_004E, S1901_C01_005E, S1901_C01_006E, 
         S1901_C01_007E, S1901_C01_008E, S1901_C01_009E, S1901_C01_010E, S1901_C01_011E) %>%
  mutate(across(-1, function(x) as.numeric(x)/100)) 

counties <- income_profile$NAME
income_profile_m <- income_profile %>% select(-NAME) %>% as.matrix()


income_similarity <- outer(1:nrow(income_profile), 1:nrow(income_profile), 
                              Vectorize(function(i, j) cosine_similarity(income_profile_m[i, ], 
                                                                         income_profile_m[j, ]))) %>%
  as.table() %>% as.data.frame() %>%
  cbind(expand_grid(county2 = counties, county1 = counties)) %>%
  mutate(county_pair = pmap_chr(list(county1, county2), 
                                ~ paste(sort(c(..1, ..2)), collapse = "+"))) %>%
  filter(county1 != county2) %>%
  select(-county1, -county2, -Var1, -Var2) %>% distinct() %>%
  separate(county_pair, into=c('county.name1', 'county.name2'), sep='\\+', remove=FALSE) %>%
  separate(county.name1, into=c('county1', 'state1'), sep=', ', remove=FALSE) %>%
  separate(county.name2, into=c('county2', 'state2'), sep=', ', remove=FALSE) 


# precipication
# https://www.ncei.noaa.gov/access/monitoring/climate-at-a-glance/county/mapping/110/pcp/202412/12/value - January-December 2022 annual precip 
precipitation <- read.csv('../data/county_precip.csv') %>%
  select(Name, State, Value) %>%
  rbind(read.csv('../data/county_precip_AK.csv') %>% select(Name, State, Value)) %>%
  rbind(read.csv('../data/county_precip_HI.csv') %>% select(Name, State, Value)) %>% 
  mutate(r_precip = (Value - min(Value))/(max(Value) - min(Value)))

# temperature
# January - December 2022 annual temperature
temperature <- read.csv('../data/county_temp.csv') %>%
  select(Name, State, Value) %>%
  rbind(read.csv('../data/county_temp_AK.csv') %>% select(Name, State, Value)) %>%
  rbind(read.csv('../data/county_temp_HI.csv') %>% select(Name, State, Value)) %>% 
  mutate(r_temp =  (Value - min(Value))/(max(Value) - min(Value)))

# population density
density <- population %>%
  left_join(county_map_data %>% select(abbr, full, county, -geom) %>% distinct(), 
            by=c('CTYNAME' = 'county', 'STNAME' = 'full')) %>%
  left_join(read.csv('../data/2024_Gaz_counties_national.txt', sep='\t') %>%
              select(NAME, USPS, ALAND_SQMI), 
            by=c('CTYNAME' = 'NAME', 'abbr' = 'USPS'))%>%
  mutate(population_density = POPESTIMATE2022/ALAND_SQMI,
         r_pd = (population_density - min(population_density))/(max(population_density)-min(population_density))) %>%
  select(CTYNAME, STNAME, r_pd)

# road density
# road length file: https://gist.github.com/erdavis1/c275afbcfffbb6bf01ade657ad0ddd65
rd <- read.csv('../data/2024_Gaz_counties_national.txt', sep='\t') %>%
  select(NAME, USPS, ALAND_SQMI) %>%
  left_join(county_map_data %>% select(abbr, full, county, -geom) %>% distinct(),
            by = c('NAME'='county', 'USPS'= 'abbr')) %>%
  left_join(read.csv('../data/road_length_by_county.csv') %>% select(State, CountyName, Suffix, Len..meters.) %>%
              group_by(State, CountyName) %>%
              summarize(length = sum(Len..meters.)/1609),
            by = c('NAME'='CountyName', 'full'='State')) %>%
  select(NAME, USPS, full, length) %>%
  left_join(read.csv('../data/2024_Gaz_counties_national.txt', sep='\t') %>%
              select(NAME, USPS, ALAND_SQMI), 
            by=c('NAME'='NAME', 'USPS'='USPS')) %>%
  mutate(road_density = length/ALAND_SQMI,
         r_rd = (road_density- min(road_density, na.rm = TRUE))/(max(road_density, na.rm = TRUE)-min(road_density, na.rm = TRUE)))

# latitude:
# dataset from https://www.census.gov/geographies/reference-files/time-series/geo/gazetteer-files.html (2024)
lat <- read.csv('../data/2024_Gaz_counties_national.txt', sep='\t') %>%
  select(NAME, USPS, INTPTLAT) %>%
  mutate(r_lat = (INTPTLAT - min(INTPTLAT))/(max(INTPTLAT) - min(INTPTLAT))) %>%
  left_join(county_map_data %>% select(abbr, full, county, -geom) %>% distinct(),
            by = c('NAME'='county', 'USPS'= 'abbr')) %>%
  select(-geom)

# life expectancy:
# source: CDC https://www.cdc.gov/nchs/data-visualization/life-expectancy/index.html
life_expectancy <- read.csv('../data/U.S._Life_Expectancy_at_Birth_by_State_and_Census_Tract_-_2010-2015.csv') %>%
  separate(County, into = c('County', 'abbr'), sep=', ') %>%
  filter(County != '(blank)',
         !is.na(abbr),
         !is.na(Life.Expectancy)) %>%
  group_by(State, County, abbr) %>%
  summarize(life_expectancy = mean(Life.Expectancy, na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(r_le = (life_expectancy - min(life_expectancy))/(max(life_expectancy) - min(life_expectancy)))


entropy_by_county <- language_diversity %>%
  select(-starts_with('Language other than'), -ends_with('"very well"')) %>%
  pivot_longer(cols = -NAME, names_to = 'category', values_to = 'value') %>%
  mutate(value = as.numeric(value)/100) %>%
  group_by(NAME) %>%
  summarise(H = sum(-value*log2(value+1e-30))) %>%
  separate(NAME, into=c('County', "State"), sep=', ', remove=FALSE) 

viz <- usmap::us_map(regions = "counties") %>%
  left_join(entropy_by_county, by=c('full'='State', 'county' = 'County')) %>%
  left_join(rd, by=c('county'='NAME', 'full', 'abbr'='USPS'))


ggplot(viz) +
  geom_sf(aes(fill=pmin(road_density, 6))) +
  scale_fill_gradient(low='yellow', high='brown')



# construct county pair dataset
# environmental factors:
# mean annual temperature, mean annual precipitation, latitude
# social factors:
# income profile, education profile, life expectancy, population density, road density

df <- income_similarity %>% rename(income_similarity = Freq) %>%
  left_join(education_similarity %>% rename(education_similarity = Freq),
            by=c("county_pair", "county.name1", "county1", "state1", "county.name2", "county2", "state2")) %>%
  left_join(temperature %>% rename(county1 = Name, state1 = State, r_temp1 = r_temp) %>% select(-Value), 
            by = c('county1', 'state1')) %>%
  left_join(temperature %>% rename(county2 = Name, state2 = State, r_temp2 = r_temp) %>% select(-Value), 
            by = c('county2', 'state2')) %>%
  left_join(precipitation %>% rename(county1 = Name, state1 = State, r_precip1 = r_precip) %>% select(-Value),
            by = c('county1', 'state1')) %>%
  left_join(precipitation %>% rename(county2 = Name, state2 = State, r_precip2 = r_precip) %>% select(-Value),
            by = c('county2', 'state2')) %>%
  left_join(lat %>% rename(county1 = NAME, state1 = full, r_lat1 = r_lat) %>% select(-INTPTLAT, -USPS),
            by = c('county1', 'state1')) %>%
  left_join(lat %>% rename(county2 = NAME, state2 = full, r_lat2 = r_lat) %>% select(-INTPTLAT, -USPS),
            by = c('county2', 'state2')) %>%
  left_join(life_expectancy %>% rename(county1 = County, state1 = State, r_le1 = r_le) %>% select(-abbr, -life_expectancy),
            by = c('county1', 'state1')) %>%
  left_join(life_expectancy %>% rename(county2 = County, state2 = State, r_le2 = r_le) %>% select(-abbr, -life_expectancy),
            by = c('county2', 'state2')) %>%
  left_join(density %>% rename(county1 = CTYNAME, state1 = STNAME, r_pd1 = r_pd),
            by = c('county1', 'state1')) %>%
  left_join(density %>% rename(county2 = CTYNAME, state2 = STNAME, r_pd2 = r_pd),
            by = c('county2', 'state2')) %>%
  left_join(rd %>% rename(county1 = NAME, state1 = full, r_rd1 = r_rd) %>% select(county1, state1, r_rd1),
            by = c('county1', 'state1')) %>%
  left_join(rd %>% rename(county2 = NAME, state2 = full, r_rd2 = r_rd) %>% select(county2, state2, r_rd2),
            by = c('county2', 'state2')) %>%
  mutate(key1 = str_extract(county1, county_regex),
         key2 = str_extract(county2, county_regex)) %>%
  left_join(timeseries_data %>% rename(key1 = Admin2, state1 = Province_State, mean_beta1 = mean_beta) %>% select(key1, state1, mean_beta1),
            by = c('key1', 'state1')) %>%
  left_join(timeseries_data %>% rename(key2 = Admin2, state2 = Province_State, mean_beta2 = mean_beta) %>% select(key2, state2, mean_beta2),
            by = c('key2', 'state2')) %>%
  left_join(entropy_by_county %>% rename(county1 = County, state1 = State, h1 = H) %>% select(county1, state1, h1),
            by = c('county1', 'state1')) %>%
  left_join(entropy_by_county %>% rename(county2 = County, state2 = State, h2 = H) %>% select(county2, state2, h2),
            by = c('county2', 'state2')) %>%
  filter(mean_beta2 > 0.06, mean_beta1 > 0.06) %>%
  mutate(
    delta_temp = r_temp2 - r_temp1,
    delta_precip = r_precip2 - r_precip1,
    delta_lat = r_lat2 - r_lat1,
    delta_income = 1 - income_similarity,
    delta_ed = 1 - education_similarity,
    delta_le = r_le2 - r_le1,
    delta_pd = r_pd2 - r_pd1,
    delta_rd = r_rd2 - r_rd1,
    delta_h = h2 - h1,
    delta_beta = mean_beta2 - mean_beta1
  ) %>% drop_na() %>%
  select(
    county_pair, 
    delta_temp, delta_precip, delta_lat, delta_income, delta_ed, delta_le,
    delta_pd, delta_rd, delta_h, delta_beta
  ) %>%
  mutate(
    across(starts_with('delta_'), ~ ifelse(delta_h < 0, -.x, .x)),
    sum_delta = abs(delta_temp) + abs(delta_precip) + abs(delta_lat) +
      abs(delta_income) + abs(delta_ed) + abs(delta_le) + abs(delta_pd) + abs(delta_rd)
    ) # make delta_h strictly greater than zero


pdf('../figs/reg_us.pdf', width=6, height=4)
ggplot(df, aes(x=delta_h, y=delta_beta)) +
  geom_hex(alpha=0.5) +
  geom_smooth(method='lm')+
  xlab(latex2exp::TeX('\\Delta\\,H')) +
  ylab(latex2exp::TeX('\\Delta\\beta')) +
  theme_bw()
dev.off()

pdf('../figs/beta_us.pdf', width=10, height=6)
ggplot(timeseries_data, aes(x=mean_beta)) +
  geom_histogram(bins=100, fill='blue') +
  theme_bw() +
  xlab(latex2exp::TeX('\\beta'))
dev.off()

lm(
  delta_beta ~ delta_h + delta_temp + delta_precip + delta_lat + delta_income + delta_ed + delta_le + delta_pd + delta_rd, 
  data = df
) %>% summary()










