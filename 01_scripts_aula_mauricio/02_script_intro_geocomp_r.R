#' ---
#' title: aula 02 - funcionamento linguagem r
#' author: Gedimar Barbosa
#' date: 2020-10-19
#' ---

#Comentários:
#Check R code style on how to write your code (standards, details, etc)

7+(7/7)+(7*7)-7

10:60 #sequencia unitaria de 10 a 60

#Atribuição de dados à objetos/ variáveis no R (<-) (Alt - é atalho para inserir <-) 

obj_10 <- 10

obj_10 #sempre confiram o objeto após a atribuição

summary(obj_10)

obj_2 <- 2

adi <- obj_10 + obj_2
adi


var <- 3*2^3
var

var1 <- 2*3^2
var1

var == var1


#Funções
#Exemplo

sum(1, 2, 3, NA)

sum(1, 2, 3, NA, na.rm = TRUE)

prod(10,2)

sum(var, var1)

#Função rep repete valores

rep(x = 1:5, times = 10) #10 vezes o conjunto todo

rep(x = 1:5, each = 10) #10 vezes o conjunto, um número de cada vez

var3 <- 100
var4 <- 300

mult <- prod(var3, var4)

mult

ln <- log(mult)

ln

#checando pacotes carregado e instalados

search()
library()

library(vegan)

library(devtools)

library(tidyverse)

#Abrir descrição detalhada de como o pacote funciona e o que pode ser feito
library(help = "vegan")

ls() #lista objetos criados

#Save workspace sempre que vc necessitar salvar os objetos criados no Global Environment
save.image("meus_objetos.RData")

#Salva só o que vc quiser
save(obj1, obj2, obj3, file = "teste.rda")

load("meus_objetos.RData")

#remover objetos do environment

rm(var) #remove o var

#remover tudo do environment
rm(list = ls())

load("meus_objetos.RData")




