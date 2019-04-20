\documentclass{article}
\usepackage[utf8]{inputenc}
\usepackage[english]{babel}


\usepackage{natbib}
\usepackage{graphicx}
\graphicspath{ {images/} }
\usepackage{mathtools}
\usepackage{listings}
\usepackage{color}
\usepackage{authblk}
\usepackage{appendix}
\usepackage{amsthm}
\usepackage{float}

\theoremstyle{definition}
\newtheorem{definition}{Definition}[section]
 
 
\definecolor{codegreen}{rgb}{0,0.6,0}
\definecolor{codegray}{rgb}{0.5,0.5,0.5}
\definecolor{codepurple}{rgb}{0.58,0,0.82}
\definecolor{backcolour}{rgb}{0.95,0.95,0.92}
\definecolor{codekeywords}{rgb}{0.337, 0.612, 0.839}
 
\definecolor{codegreen}{rgb}{0,0.6,0}
\definecolor{codegray}{rgb}{0.5,0.5,0.5}
\definecolor{codepurple}{rgb}{0.58,0,0.82}
\definecolor{backcolour}{rgb}{0.95,0.95,0.92}
 
\lstdefinestyle{mystyle}{
    backgroundcolor=\color{backcolour},   
    commentstyle=\color{codegreen},
    keywordstyle=\color{magenta},
    numberstyle=\tiny\color{codegray},
    stringstyle=\color{codepurple},
    basicstyle=\footnotesize,
    breakatwhitespace=false,         
    breaklines=true,                 
    captionpos=b,                    
    keepspaces=true,                 
    numbers=left,                    
    numbersep=5pt,                  
    showspaces=false,                
    showstringspaces=false,
    showtabs=false,                  
    tabsize=2
}
 
\lstset{style=mystyle}


\newcommand\inner[2]{\langle #1, #2 \rangle}
\newcommand\norm[1]{\| #1 \|}


\title{Writing a Basic Ray Tracer}
\author{Vassillen Chizhov}
\affil{Saarland University, IMPRS}
\date{April 2019}


<p align="center"><img src="/tex/46e7fcb2cb755d1452a04c0d062815ab.svg?invert_in_darkmode&sanitize=true" align=middle width=1790.9885399999998pt height=26677.47930045pt/></p>

