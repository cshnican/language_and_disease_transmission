library(tidyverse)

hdi <- read.csv('../data/GDL-Subnational-HDI-data.csv') %>%
  pivot_longer(starts_with('X'), names_to = 'year', values_to = 'HDI') %>%
  mutate(year = as.numeric(str_sub(year, 2, -1)))

temperature <- read.csv('../data/GDL-Yearly-Average-Surface-Temperature-(ºC)-data.csv') %>%
  pivot_longer(starts_with('X'), names_to = 'year', values_to = 'mean_temperature') %>%
  mutate(year = as.numeric(str_sub(year, 2, -1)))

precip <- read.csv('../data/GDL-Total-Yearly-Precipitation-(m)-data.csv') %>%
  pivot_longer(starts_with('X'), names_to = 'year', values_to = 'mean_precipitation') %>%
  mutate(year = as.numeric(str_sub(year, 2, -1)))

elevation <- read.csv('../data/elevations.csv', sep='\t')

language_diversity <- read.csv('../data/ethonologue_diversity.csv', sep='\t') %>%
  mutate(diversity_index = as.numeric(diversity_index))

# from Worldbank
population <- read.csv('../data/populations.csv') %>%
  pivot_longer(starts_with('X'), names_to = 'year', values_to = 'population') %>%
  mutate(year = as.numeric(str_sub(year, 2, -1))) %>%
  filter(year == 2022) %>%
  dplyr::select(Country.Name, population)

# population density, source: https://data.worldbank.org/indicator/EN.POP.DNST
density <- read.csv('../data/API_EN.POP.DNST_DS2_en_csv_v2_89.csv') %>% 
  pivot_longer(starts_with('X'), names_to = 'year', values_to = 'population_density') %>%
  mutate(year = as.numeric(str_sub(year, 2, -1))) %>%
  filter(year == 2022) %>%
  rename(Country = Country.Name) %>%
  dplyr::select(Country, population_density)


# source: https://prosperitydata360.worldbank.org/en/indicator/WEF+TTDI+ROADDENS
roaddensity <- read.csv('../data/WEF-TTDI.csv') %>%
  filter(Indicator.ID == 'WEF.TTDI.ROADDENS', Attribute.1 == 'Value') %>%
  rename(road_density = X2021,
         Country = Economy.Name) %>%
  select(Country, road_density)

# source: Worldbank API
geo <- worldbank::wb_country() %>% 
  rename(Country = country_name) %>%
  dplyr::select(Country, latitude) %>%
  filter(!is.na(latitude))
  


# https://github.com/CSSEGISandData/COVID-19/blob/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv
# Methods from Kirkeby et al. (2017)
timeseries_data <- read.csv('../data/time_series_covid19_confirmed_global.csv') %>% 
  pivot_longer(starts_with('X'), names_to = 'date', values_to = 'cumulative_cases') %>%
  mutate(date = str_sub(date, 2, -1),
         date = as.Date(date, '%m.%d.%y')) %>%
  separate(date, into=c('year', 'month', 'day'), remove=FALSE) %>%
  group_by(Country.Region, date, month, day, year) %>%
  summarize(cumulative_cases = sum(cumulative_cases)) %>%
  left_join(population, by=c('Country.Region' = 'Country.Name')) %>% ungroup() %>%
  group_by(Country.Region) %>%
  mutate(IN = c(0, diff(cumulative_cases)),
         S = population - cumulative_cases,
         beta = (-log(1 - IN/(cumulative_cases + 1e-30) + 1e-30))/(1*S/population)) %>%
  summarize(mean_beta = mean(beta, na.rm = TRUE))
  


pdf('../figs/transmission_rate.pdf', width=8, height=8)
ggplot(timeseries_data %>% filter(!is.nan(mean_beta)) %>%
         mutate(half=mean_beta < 0.07125), 
       aes(x=reorder(Country.Region, mean_beta), y=mean_beta, fill=mean_beta)) +
  facet_wrap(~half, scales = 'free_y') +
  geom_col()+
  coord_flip() +
  xlab('Country') +
  theme_classic(8) +
  theme(strip.text.x = element_blank())
dev.off()

df <- hdi %>%
  filter(Level == 'National') %>%
  left_join(temperature, by = c('Country', 'Continent', 'ISO_Code', 'Level', 'GDLCODE', 'Region', 'year')) %>%
  left_join(precip %>% distinct() %>% filter(!is.na(mean_precipitation)), 
            by = c('Country', 'Continent', 'ISO_Code', 'Level', 'GDLCODE', 'Region', 'year')) %>%
  left_join(elevation, by = 'Country') %>%
  left_join(language_diversity %>% dplyr::select(Country, diversity_index), by = 'Country') %>%
  left_join(roaddensity, by='Country') %>%
  left_join(geo, by = 'Country') %>%
  left_join(density, by = 'Country') %>%
  filter(!is.na(mean_temperature),
         !is.na(HDI),
         !is.na(Elevation),
         !is.na(diversity_index),
         !is.na(road_density),
         !is.na(latitude),
         !is.na(population_density))%>%
  mutate(r_HDI = scale(HDI)[,1],
         r_temp = scale(mean_temperature)[,1],
         r_precip = scale(mean_precipitation)[,1],
         r_elevation = scale(Elevation)[,1],
         r_road_density = scale(road_density)[,1],
         r_latitude = scale(abs(latitude))[,1],
         r_population_density = scale(population_density)[,1])
  


df_2022 <- df %>% filter(year == 2022)

# calculate pairwise difference in HDI on a national level
GDLCODE1 <- df_2022$GDLCODE[df_2022$Level == 'National'] %>% unique()
GDLCODE2 <- GDLCODE1

# exclude countries with "abnormal" beta (i.e. China, North Korea, Thailand, USA, Japan)
to_exclude <- c('China', "Korea, Dem. People's Rep.", "United States", "Thailand", "Japan")

delta_national <- expand_grid(GDLCODE1, GDLCODE2) %>%
  filter(GDLCODE1 != GDLCODE2) %>%
  left_join(df_2022 %>%
              dplyr::select(GDLCODE, Country, Region, 
                            r_HDI, 
                            r_temp, r_precip, r_elevation, r_latitude,
                            r_road_density, r_population_density,
                            diversity_index), 
            by=c('GDLCODE1' = 'GDLCODE')) %>%
  rename(HDI1=r_HDI, 
         temp1 = r_temp, precip1 = r_precip, elevation1 = r_elevation, latitude1 = r_latitude,
         road_density1 = r_road_density, population_density1 = r_population_density,
         country1 = Country, 
         region1 = Region, diversity_index1 = diversity_index) %>%
  left_join(df_2022 %>% 
              dplyr::select(GDLCODE, Country, Region, r_HDI, 
                            r_temp, r_precip, r_elevation, r_latitude,
                            r_road_density, r_population_density,
                            diversity_index), 
            by=c('GDLCODE2' = 'GDLCODE')) %>%
  rename(HDI2=r_HDI, 
         temp2= r_temp, precip2 = r_precip, elevation2 = r_elevation, latitude2 = r_latitude,
         road_density2 = r_road_density, population_density2 = r_population_density,
         country2 = Country, 
         region2 = Region, elevation2 = r_elevation, diversity_index2 = diversity_index) %>%
  left_join(timeseries_data %>% filter(!Country.Region %in% to_exclude), by=c('country1'='Country.Region')) %>%
  rename(mean_beta1 = mean_beta) %>%
  left_join(timeseries_data  %>% filter(!Country.Region %in% to_exclude), by=c('country2'='Country.Region')) %>%
  rename(mean_beta2 = mean_beta) %>%
  mutate(delta_HDI = HDI2 - HDI1,
         delta_temp = temp2 - temp1,
         delta_precip = precip2 - precip1,
         delta_elevation = elevation2 - elevation1,
         delta_latitude = latitude2 - latitude1,
         delta_rd = road_density2 - road_density1,
         delta_pd = population_density2 - population_density1,
         delta_diversity_index = diversity_index2 - diversity_index1,
         delta_beta = mean_beta2 - mean_beta1) %>%
  filter(!is.na(delta_beta)) %>%
  dplyr::select(country1, region1, GDLCODE1, country2, region2, GDLCODE2, 
                delta_HDI, 
                delta_temp, delta_precip, delta_elevation, delta_latitude,
                delta_rd, delta_pd,
                delta_diversity_index, 
                delta_beta) %>%
  mutate(sum_delta = abs(delta_HDI) +
          abs(delta_temp) + abs(delta_precip) + abs(delta_elevation) + abs(delta_latitude) +
           abs(delta_rd) + abs(delta_pd)) %>%
  filter(country1 != country2,
         delta_diversity_index >=0) %>%
  mutate(country_pair = pmap_chr(list(country1, country2), 
                                 ~ paste(sort(c(..1, ..2)), collapse = "+"))) %>%
  group_by(country_pair) %>%
  filter(row_number() == 1) %>% ungroup()


final <- delta_national

pdf('../figs/diversity_index.pdf', width=6, height=12)
ggplot(final %>%
         rename(`Sum of diff in metrics` = sum_delta),
       aes(x = reorder(country_pair, delta_diversity_index), 
           y = delta_diversity_index, fill = `Sum of diff in metrics`)) +
  geom_col() +
  coord_flip() +
  xlab('Country pair with the most similar conditions') +
  ylab(latex2exp::TeX('\\Delta LDI'))
dev.off()



pdf('../figs/reg.pdf', width=6, height=4)
ggplot(final,#delta_national %>% filter(abs(delta_beta) < 0.04),
       aes(x = delta_diversity_index, y = delta_beta)) +
  geom_hex()+
  #ggrepel::geom_text_repel(aes(label = country_pair)) +
  geom_smooth(method='lm')+
  xlab(latex2exp::TeX('\\Delta\\,LDI')) +
  ylab(latex2exp::TeX('\\Delta\\beta')) +
  theme_bw()
dev.off()


lm(
  delta_beta ~ delta_diversity_index + delta_pd + delta_rd + delta_elevation + delta_HDI + delta_precip + delta_temp + delta_latitude, 
  data = final) %>% summary()



# # without the delta
# national <- df_2022 %>%
#   left_join(timeseries_data %>% filter(!Country.Region %in% to_exclude), by=c('Country'='Country.Region')) %>%
#   filter(!is.na(mean_beta))
# 
# ggplot(national, aes(x=diversity_index, y=mean_beta)) +
#   geom_point() +
#   geom_smooth(method='lm')
# 
# lm(
#   mean_beta ~ diversity_index + 
#     r_temp + r_precip + r_elevation + r_latitude +
#     r_HDI + r_population_density + r_road_density,
#   data = national) %>% summary()


