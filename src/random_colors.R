scPalette <- function(n) {
  colorSpace <- c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#F29403','#F781BF','#BC9DCC','#A65628','#54B0E4','#222F75','#1B9E77','#B2DF8A',
                           '#E3BE00','#FB9A99','#E7298A','#910241','#00CDD1','#A6CEE3','#CE1261',"#69AED4","#B7B44B","#D96B7C","#69BE54","#939325",
                                           "#E5E576","#127845","#761A55","#6A32A0","#77471F",
                                           "#292774","#AA7947",
                                           "#DFAC7B","#7DC783","#0A8CED","#FEEE75","#E3AD28",
                                           "#5576A6",
                                           "#61DDF2","#4888C7","#71ACD8","#BD80B8",
                                           "#34A07C","#AC4689",
                                           "#7472B5",
                                            '#FBFAD2','#CBE5C5',"#E83D36", "#62B8BF","#B0CFE0",
                                            "#B6D88F","#FFDCAC","#DBDCDA","#F9DDE9", #Inh VIP
                                            '#2FA147','#8466AF',"#EFE3CD", # Inh LAMP5
                                            "#F8B6B4","#90D2D7","#B85E2C","#3581BA", # Inh SST
                                            "#F5999D","#A27833","#FBBE6F","#BD96EA",# Inh PVALB
                                            "#DFCCE4", '#FDD3BA',# Inh PAX6
                                            "#FFDEDE","#FFB0BD","#D7517C","#991343", #Ast
                                            "#F9CFA5","#f79c5e","#b55a00", # Oli
                                            "#FFE2CF","#b38766","#75421f",# OPC
                                            "#E1CAE5","#573B88", #Mic
                                            "#7F0009","#6D0775", #Mac, T
                                            "#E283B5","#F2C0D5","#6EC5A6","#A45627","#F38C64","#8BA0CC",
                  '#5E4FA2','#8CA77B','#00441B','#DEDC00','#B3DE69','#8DD3C7','#999999')
                           if (n <= length(colorSpace)) {
                             colors <- colorSpace[1:n]
                           } else {
                             colors <- grDevices::colorRampPalette(colorSpace)(n)
                           }
  return(colors)
}
