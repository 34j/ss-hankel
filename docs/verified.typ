#import "@preview/physica:0.9.2": *
#import "@preview/cetz:0.3.0"


// simple page setup
#let mainfont = ("BIZ UDPMincho")
#set text(font: mainfont, lang: "ja", size: 9pt)
#set page(numbering: "1")
#set math.equation(numbering: "1.")
#let zeros = $op("Zeros")$
#let poles = $op("Poles")$

$TT$: 単位円周

#definition[

  $
  mu_k := 1/(2 pi i) integral_TT
  $
]
