;; -*- lexical-binding: t; -*-

(TeX-add-style-hook
 "presentation"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("beamer" "14pt" "t")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("pgfplots" "") ("fontenc" "T1") ("graphicx" "") ("amsmath" "") ("amssymb" "") ("gensymb" "") ("babel" "slovene") ("biblatex" "sorting=none" "style=numeric" "sortcites=true") ("siunitx" "")))
   (add-to-list 'LaTeX-verbatim-environments-local "semiverbatim")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "href")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "beamer"
    "beamer10"
    "fontenc"
    "graphicx"
    "amsmath"
    "amssymb"
    "gensymb"
    "babel"
    "biblatex"
    "siunitx")
   (LaTeX-add-bibliographies
    "../source/refs"))
 :latex)

