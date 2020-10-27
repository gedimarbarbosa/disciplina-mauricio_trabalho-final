#' ---
#' title: aula 04 introducao ao tidyverse
#' author: mauricio vancine
#' date: 2020-10-20
#' ---


# packages ----------------------------------------------------------------
library(tidyverse)
library(here)
library(readxl)
library(writexl)
library(lubridate)

# topics ------------------------------------------------------------------
# 4.1 tidyverse
# 4.2 magrittr (pipe - %>%)
# 4.3 readr
# 4.4 readxl e writexl
# 4.5 tibble
# 4.6 tidyr
# 4.7 dplyr
# 4.8 stringr
# 4.9 forcats
# 4.10 lubridate
# 4.11 purrr

# 4.1 tidyverse -----------------------------------------------------------
# instalar o pacote
# install.packages("tidyverse")

# carregar o pacote
library(tidyverse)

# list all packages in the tidyverse 
tidyverse::tidyverse_packages(include_self = TRUE)

# 4.2 magrittr (pipe - %>%) -----------------------------------------------
# sem pipe
sqrt(sum(1:100))

# com pipe (funciona como se fosse uma cascata de operações) -- Atalho ctrl + shift + m
1:100 %>% 
  sum() %>% 
  sqrt()

# fixar amostragem
set.seed(42)

# sem pipe
ve <- sum(sqrt(sort(log10(rpois(100, 10)))))
ve

# fixar amostragem
set.seed(42)

# com pipe
ve <- rpois(100, 10) %>% 
  log10() %>%
  sort() %>% 
  sqrt() %>% 
  sum()
ve  

# exercicio 09 ------------------------------------------------------------

log10(cumsum(1:100))
sum(sqrt(rnorm(100)))
prod(sort(sample(1:10, 10000, rep = TRUE)))

1:100 %>% 
  cumsum %>% 
  log10

rnorm(100) %>% 
  abs() %>% 
  sqrt() %>% 
  sum()

1:10 %>% 
sample(10000, rep = TRUE) %>% 
  sort %>% 
  sum #não precisa do ()

# 4.3 readr ---------------------------------------------------------------
# install
# install.packages("here")

# load
library(here)

# confer #procura onde está o project
here::here()

# create a .here file
# here::set_here() #cria um here no diretório

here::set_here()

# formato .csv
# import sites
si <- readr::read_csv(here::here("03_dados", "tabelas", "ATLANTIC_AMPHIBIANS_sites.csv"))
si

# import sites #caminho relativo sem usar o here
si <- readr::read_csv("./03_dados/tabelas/ATLANTIC_AMPHIBIANS_sites.csv")
si

# formato .txt
# import sites
si <- readr::read_tsv(here::here("03_dados", "tabelas", "ATLANTIC_AMPHIBIANS_sites.txt"))
si

# 4.4 readxl e writexl ----------------------------------------------------
# import .xlsx
# install.packages("readxl")
library("readxl")

# export .xlsx
# install.packages("writexl")
library("writexl")

# import sites
si <- readxl::read_xlsx(here::here("03_dados", "tabelas", "ATLANTIC_AMPHIBIANS_sites.xlsx"), 
                        sheet = 1)
si

# import sites
si <- readr::read_csv(here::here("03_dados", "tabelas", "ATLANTIC_AMPHIBIANS_sites.csv"))
si

# import species
sp <- readr::read_csv(here::here("03_dados", "tabelas", "ATLANTIC_AMPHIBIANS_species.csv"))
sp

# 4.5 tibble --------------------------------------------------------------
# view the sites data
glimpse(si)

# view the species data
glimpse(sp)

# tibble vs data.frame
# 1. nunca converte um tipo character como factor - 
df <- data.frame(ch = c("a", "b"), nu = 1:2)
str(df)

tb <- tibble(ch = c("a", "b"), nu = 1:2)
glimpse(tb)

# 2. a indexacao com colchetes sempre retorna um tibble
df_ch <- df[, 1]
class(df_ch)

tb_ch <- tb[, 1]
class(tb_ch)

# indexacao pelo nome devolve um vetor
tb_ch <- tb$ch
class(tb_ch)

# 3. nao faz correspondencia parcial, retorna NULL se a coluna nao existe com o nome especificado
df$c 
tb$c

# 4.6 tidyr ---------------------------------------------------------------
# funcoes
# 1 unite(): junta dados de multiplas colunas em uma
# 2 separate(): separa caracteres em mulplica colunas
# 3 drop_na(): retira linhas com NA
# 4 replace_na(): substitui NA
# 5 spread(): long para wide
# 6 gather(): wide para long
  
# 1 unite
# sem pipes
si_unite <- tidyr::unite(si, "lat_lon", latitude:longitude, sep = ",")
si_unite$lat_lon

# com pipes
si_unite <- si %>% 
  tidyr::unite("lat_lon", latitude:longitude, sep = ",")
si_unite$lat_lon
  
# 2 separate
# separar os dados de "period" em quatro colunas dos seus valores
si_separate <- si %>% 
  tidyr::separate("period", c("mo", "da", "tw", "ni"), remove = FALSE) #se colocarmos TRUE ela remove a coluna period
si_separate[, c(1, 9:13)]

# 3 separate_rows()
# Separar os dados de "period" na mesma coluna e repetindo os valores
si_separate_row <- si %>% 
  tidyr::separate_rows("period")
si_separate_row[, c(1, 9:13)]

si_separate_row <- si %>% 
  tidyr::separate_rows("period")
si_separate_row[, c(1, 9:13)]

# 4 drop_na()
# remove as linhas com NA de todas as colunas
si_drop_na <- si %>% 
  tidyr::drop_na()
si_drop_na

# remove as linhas com NA da coluna "active_methods"
si_drop_na <- si %>% 
  tidyr::drop_na(active_methods)
si_drop_na

# 5 replace_na()
# Substituir os NAs da coluna "active_methods" por 0 
si_replace_na <- si %>% 
  tidyr::replace_na(list(active_methods = 0))
si_replace_na

# 6 pivot_wider():  long para wide ######ISSO É IMPORTANTE

# 1. id_cols: variavel id
# 2. names_from: variavel categorica que ira definir os nomes das colunas
# 3. values_from: variavel numerica que ira preencher os dados
# 4. values_fill: valor para preencher os NAs

# sites
si[, c("id", "state_abbreviation", "species_number")]

si_wide <- si %>% 
  tidyr::pivot_wider(id_cols = id, 
                     names_from = state_abbreviation,
                     values_from = species_number, 
                     values_fill = list(species_number = 0))
si_wide

# species
sp[1:1000, c("id", "species", "individuals")]

sp_wide <- sp[1:1000, ] %>% 
  tidyr::replace_na(list(individuals = 0)) %>% 
  tidyr::pivot_wider(id_cols = id, 
                     names_from = species, 
                     values_from = individuals, 
                     values_fill = list(individuals = 0))
sp_wide

# 7 pivot_longer(): wide para long

# 1. cols: coluna do id
# 2. names_to: nome da coluna que receberá os nomes
# 3. values_to: nome da coluna que receberá os valores

sP_long <- sP_wide %>% 
  tidyr::pivot_longer(cols = -id, 
                      names_to = "SPECIES", 
                      values_to = "species_number")
sP_long

sp_long <- sp_wide %>% 
  tidyr::pivot_longer(cols = -id, 
                      names_to = "species", 
                      values_to = "individuals") %>% 
  dplyr::filter(individuals != 0)
sp_long

# exercicio 10 ------------------------------------------------------------

si_un <- si %>% 
  tidyr::unite("local_total", country:site, sep = ", ")
si_un$local_total

# exercicio 11 ------------------------------------------------------------

family_wide <- sp[1:1000, ] %>% 
  tidyr::replace_na(list(individuals = 0)) %>% 
  tidyr::drop_na(family) %>% 
  tidyr::pivot_wider(id_cols = id, 
                     names_from = family, 
                     values_from = individuals, 
                     values_fn = list(individuals = sum),
                     values_fill = list(individuals = 0))
family_wide

# 4.7 dplyr ---------------------------------------------------------------
# funcoes
# 1 select(): seleciona colunas pelo nome gerando um tibble
# 2 pull(): seleciona uma coluna como vetor
# 3 rename(): muda o nome das colunas
# 4 mutate(): adiciona novas colunas ou adiciona resultados em colunas existentes
# 5 arrange(): reordenar as linhas com base nos valores de colunas
# 6 filter(): seleciona linhas com base em valores
# 7 distinct(): remove linhas com valores repetidos com base nos valores de colunas
# 8 slice(): seleciona linhas pelos numeros
# 9 n_sample(): amostragem aleatoria de linhas
# 10 summarise(): agrega ou resume os dados atraves de funcoes, podendo considerar valores das colunas
# 11 *_join(): junta dados de duas tabelas atraves de uma coluna chave


# 1 select
# seleciona colunas pelo nome
si_select <- si %>% 
  dplyr::select(id, longitude, latitude)
si_select

# nao seleciona colunas pelo nome
si_select <- si %>% 
  select(-c(id, longitude, latitude))
si_select

#  starts_with(), ends_with() e contains()
si_select <- si %>% 
  select(contains("sp"))
si_select

# 2 pull
# coluna para vetor
si_pull <- si %>% 
  pull(id)
si_pull

si_pull <- si %>% 
  pull(species_number)
si_pull

# 3 rename
si_rename <- si %>%
  rename(sp = species_number)
si_rename

# 4 mutate
si_mutate <- si %>% 
  mutate(record_factor = as.factor(record))
si_mutate$record_factor

si_mutate

# 5 arrange
si_arrange <- si %>% 
  arrange(species_number)
si_arrange

si_arrange <- si %>% 
  arrange(desc(species_number))
si_arrange

si_arrange <- si %>% 
  arrange(-species_number)
si_arrange

# 6 filter
si_filter <- si %>% 
  filter(species_number > 5)
si_filter

si_filter <- si %>% 
  filter(between(species_number, 1, 5))
si_filter

si_filter <- si %>% 
  filter(is.na(passive_methods))
si_filter

si_filter <- si %>% 
  dplyr::filter(!is.na(passive_methods)) #! significa uma negação
si_filter

si_filter <- si %>% 
  filter(!is.na(active_methods) & !is.na(passive_methods))
si_filter

si_filter <- si %>% 
  filter(species_number > 5 & state_abbreviation == "BR-SP") #& E
si_filter

si_filter <- si %>% 
  filter(species_number > 5 | state_abbreviation == "BR-SP")#| OU 
si_filter

# 7 distinct
si_distinct <- si %>% 
  distinct(species_number)
si_distinct

si_distinct <- si %>% 
  distinct(species_number, .keep_all = TRUE)
si_distinct

# 8 slice 
si_slice <- si %>% 
  slice(1:10)
si_slice

# 9 n_sample 
si_sample_n <- si %>% 
  sample_n(100)
si_sample_n

# 10 summarise
si_summarise <- si %>% 
  summarise(mean_sp = mean(species_number), 
            sd_sp = sd(species_number))
si_summarise

si_summarise_group <- si %>% 
  group_by(country) %>% #agrupando primeiro por país
  summarise(mean_sp = mean(species_number), 
            sd_sp = sd(species_number))
si_summarise_group

# 11 join #no meu trabalho, posso usar full join da para juntar as matrizes (mantem tudo)
si_coord <- si %>% 
  select(id, longitude, latitude)
si_coord 

sp_join <- sp %>% 
  left_join(si_coord, by = "id")
sp_join

colnames(sp_join)

sp_join %>% 
  dplyr::select(species, longitude, latitude)

# sufixos _at(), _if(), _all(): realiza operações dependente de confições
sp_wide_rename <- sp_wide %>% 
  dplyr::rename_at(vars(contains(" ")), 
                   list(~stringr::str_replace_all(., " ", "_"))) %>% 
  dplyr::rename_all(list(~stringr::str_to_lower(.)))
sp_wide_rename

#  manipular os dados de forma mais simples
da <- si %>% 
  select(id, state_abbreviation, species_number)
da

da <- si %>% 
  select(id, state_abbreviation, species_number) %>% 
  filter(species_number > 5)
da

da <- si %>% 
  select(id, state_abbreviation, species_number) %>% 
  filter(species_number > 5) %>% 
  group_by(state_abbreviation) %>% 
  summarise(nsp_state_abb = n())
da

da <- si %>% 
  select(id, state_abbreviation, species_number) %>% 
  filter(species_number > 5) %>% 
  group_by(state_abbreviation) %>% 
  summarise(nsp_state_abb = n()) %>% 
  arrange(nsp_state_abb)
da


# exercicio 12 ------------------------------------------------------------

si %>% 
  dplyr::mutate(alt_log = log10(altitude), 
                tem_log = log10(temperature),
                pre_log = log10(precipitation))


# exercicio 13 ------------------------------------------------------------

si %>% 
  dplyr::arrange(-altitude)

# exercicio 14 ------------------------------------------------------------

si %>% 
  dplyr::filter(altitude > 1000,
                temperature < 15,
                precipitation > 1000)

# exercicio 15 ------------------------------------------------------------

si %>% 
  dplyr::filter(species_number>15) %>% 
  dplyr::sample_n(200)

# 4.8 stringr -------------------------------------------------------------
# comprimento
stringr::str_length("abc")

# extrai
stringr::str_sub("abc", 3)

# inserir espaco em branco
stringr::str_pad("abc", width = 4, side = "left")
stringr::str_pad("abc", width = 4, side = "right")

# remover espaco em branco do comeco, final ou ambos
stringr::str_trim(" abc ")

# minusculas e maiusculas
stringr::str_to_upper("abc")
stringr::str_to_title("abc")
stringr::str_to_title("aBc")

# ordenacao
le <- sample(letters, 26, rep = TRUE)
le

stringr::str_sort(le)
stringr::str_sort(le, dec = TRUE)

# extrair
stringr::str_extract("abc", "b")

# substituir
stringr::str_replace("abc", "a", "y")

# separacao
stringr::str_split("a-b-c", "-")

# 4.9 forcats -------------------------------------------------------------
# fixar amostragem
set.seed(42)

# as_factor(): cria fatores
fa <- sample(c("alto", "medio", "baixo"), 30, rep = TRUE) %>% 
  forcats::as_factor()
fa

# fct_recode(): muda o nome dos níveis
fa_recode <- fa %>% 
  forcats::fct_recode(a = "alto", m = "medio", b = "baixo")
fa_recode

# fct_rev(): inverte os niveis
fa_rev <- fa_recode %>% 
  forcats::fct_rev()
fa_rev

# fct_relevel(): especifica a classificacao de um nivel
fa_relevel <- fa_recode %>% 
  forcats::fct_relevel(c("a", "m", "b"))
fa_relevel

# fct_inorder(): ordem em que aparece
fa_inorder <- fa_recode %>% 
  forcats::fct_inorder()
fa_inorder

# fct_infreq: ordem (decrescente) de frequencia
fa_infreq <- fa_recode %>% 
  forcats::fct_infreq()
fa_infreq

# fct_lump(): agregacao de niveis raros em um nivel
fa_lump <- fa_recode %>% 
  forcats::fct_lump()
fa_lump


# 4.10 lubridate ----------------------------------------------------------
#Bom para trabalhar com questões temporais, envolvendo data, hora, etc

# install
install.packages("lubridate")

# load
library(lubridate)

# string
data_string <- "2020-04-24"
data_string
class(data_string)

# criar um objeto com a classe data
data_date <- lubridate::date(data_string)
data_date
class(data_date)

# criar um objeto com a classe data
data_date <- lubridate::as_date(data_string)
data_date
class(data_date)

# string
data_string <- "20-10-2020"
data_string
class(data_string)

# criar um objeto com a classe data - dia/ mes/ ano
data_date <- lubridate::dmy(data_string)
data_date
class(data_date)

# formatos
lubridate::dmy(20102020)
lubridate::dmy("20102020")
lubridate::dmy("20/10/2020")
lubridate::dmy("20.10.2020")

# especificar horarios
lubridate::dmy_h(2010202020)
lubridate::dmy_hm(201020202035)
lubridate::dmy_hms(20102020203535)

# criar
data <- lubridate::dmy_hms(20102020203535)
data

# extrair
lubridate::second(data)
lubridate::day(data)
lubridate::month(data)
lubridate::wday(data) #dia da semana
lubridate::wday(data, label = TRUE) #dia da semana

# inlcuir
lubridate::hour(data) <- 13
data

# extrair a data no instante da execucao
lubridate::today() 

# extrair a data e horario no instante da execucao
lubridate::now()

# agora
agora <- lubridate::ymd_hms(lubridate::now(), tz = "America/Sao_Paulo")
agora

# verificar os tz
OlsonNames()

# que horas sao em...
lubridate::with_tz(agora, tzone = "GMT")
lubridate::with_tz(agora, tzone = "Europe/Stockholm")  

# altera o fuso sem mudar a hora
lubridate::force_tz(agora, tzone = "GMT")

# datas
inicio_r <- lubridate::dmy("30-11-2011")
inicio_r

hoje_r <- lubridate::today()
hoje_r

# intervalo
r_interval <- lubridate::interval(inicio_r, hoje_r)
r_interval
class(r_interval)

# outra forma de definir um intervalo: o operador %--%
r_interval <- lubridate::dmy("30-11-2011") %--% lubridate::today() 
namoro_interval <- lubridate::dmy("25-06-2008") %--% lubridate::today()   

# verificar sobreposicao
lubridate::int_overlaps(r_interval, namoro_interval)

# somando datas
inicio_r + lubridate::ddays(1)
inicio_r + lubridate::dyears(1)

# criando datas recorrentes
reunioes <- lubridate::today() + lubridate::weeks(0:10)
reunioes

# duracao de um intervalo 
r_interval <- inicio_r %--% lubridate::today()
r_interval

# transformacoes
r_interval / lubridate::dyears(1)
r_interval / lubridate::ddays(1)

# total do periodo estudando r
lubridate::as.period(r_interval)

# tempo de namoro
lubridate::as.period(namoro_interval)

# 4.11 purrr --------------------------------------------------------------
#Pacote importante para meus trabalhos

# lista
x <- list(1:5, c(4, 5, 7), c(1, 1, 1), c(2, 2, 2, 2, 2))
x

purrr::map(x, sum)
purrr::map_dbl(x, sum)
purrr::map_chr(x, paste, collapse = " ")

x <- list(3, 5, 0, 1)
y <- list(3, 5, 0, 1)

purrr::map2_dbl(x, y, prod)

x <- list(3, 5, 0, 1)
y <- list(3, 5, 0, 1)
z <- list(3, 5, 0, 1)

purrr::pmap_dbl(list(x, y, z), prod)

mean_var <- si %>% 
  dplyr::select(species_number, altitude) %>% 
  purrr::map_dbl(mean)
mean_var

sd_var <- si %>% 
  dplyr::select(species_number, altitude) %>% 
  purrr::map_dbl(sd)
sd_var

# end ---------------------------------------------------------------------