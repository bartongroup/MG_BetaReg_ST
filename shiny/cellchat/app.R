library(qs2)
library(shiny)
library(CellChat)

obj <- qs2::qs_read("../data/cellchat_1.qs")

CellChat::runCellChatApp(obj)