;; -*- lexical-binding: t; -*-

(TeX-add-style-hook
 "presentation"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("beamer" "14pt" "t")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("fontenc" "T1") ("graphicx" "") ("amsmath" "") ("amssymb" "") ("gensymb" "") ("babel" "slovene") ("siunitx" "") ("pgfplots" "")))
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
    "siunitx"
    "pgfplots"))
 :latex)

