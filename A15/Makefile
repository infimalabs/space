#!/usr/bin/make -f

# Group/Version.
override g := A15
override v := $(or $v,draft)

# Figures.
pngs := $(basename $(notdir $(wildcard figure/*.png.txt)))

# Constants.
gitmod := $(shell git log -1 --pretty=format:%ct)
giturl := $(subst git@,https://,$(subst :,/,$(subst .git,,$(shell git config remote.origin.url))))

# Group sources.
toolpy = $g.py

# Build dirs.
bldall = build
bldone = $(bldall)/$v
bldtmp = $(bldone)/tmp

# Build products.
bldpdf = $(bldone)/$g.pdf
blddocs = $(bldone)/$g.yml
bldpdfs = $(bldtmp)/$g.pdf
bldpngs = $(addprefix $(bldone)/,$(pngs))

# Build tooling.
pdfcpu = pdfcpu
pandoc = pandoc
python = python3
mkdir = mkdir -v
open = open
cp = cp -v
mv = mv -v
rm = rm -v

# Version commands.
finalstampc = :

# Default commands.
figc = SAVEFIG=$@ $(python) $(addprefix ./,$^)
stampv1 = sc: 1 abs, points: 42, ro: 0, op: .75, aligntext: r, position: bl, offset: 1.5625 1, color: \#040404, \
	strokecolor: \#1a1a1a, mode: 2, margin: 4, backgroundcolor: \#ffb400, border: 1 round \#ffb400
stampv2 = sc: .25 abs, points: 42, ro: 0, op: .75, aligntext: l, position: bl, offset: 2.875 1.125, color: \#040404, margins: 2
stampc = $(or $($(v)stampc),$(pdfcpu) watermark add -mode text -u in -- '$(1)' '$(2)' '$(3)')
latexv = \def\setgraphicspath{\graphicspath{$(addprefix {,$(addsuffix },$(sort $(^D))))}}
latexc = latexmk -outdir=$(bldtmp) -jobname=$(@F).tmp -usepretex='$(latexv)' -pdflua -quiet -halt-on-error -Werror $< </dev/null\
				 && $(mv) $@.tmp.pdf $@
pdfc = $(pdfcpu) merge -m create $@.tmp.pdf $^ \
			 && $(call stampc,$v,$(stampv1),$@.tmp.pdf) \
			 && $(call stampc,$(giturl),$(stampv2),$@.tmp.pdf) \
			 && $(pdfcpu) optimize $@.tmp.pdf $@
docsc = $(pandoc) $g.tex --fail-if-warnings --standalone --wrap=none --template=$(@F) --output=$@ --default-image-extension=png \
			 --resource-path=$(bldtmp):$(bldone) --to=html --mathml --citeproc --toc \
			 -M toc-title:Table\ of\ Contents -M date:$(gitmod)

# Ensure dirs.
$(bldone) $(bldtmp): ; $(mkdir) -p $@

# Compile figures.
$(bldone)/%.png: $(toolpy) figure/%.png.txt | $(bldone); $(figc)

# Compile content.
$(bldtmp)/%.pdf: %.tex $(bldpngs) | $(bldtmp); $(latexc)

# Compile website.
$(blddocs): $g.tex $g.yml $(bldpngs) | $(bldpdf); $(docsc)

# Compile document.
$(bldpdf): $(bldpdfs); $(pdfc)

# Group targets.
build: $(blddocs) $(bldpdf)

# Group tests.
t test check: build; @echo 'TODO($<): $@'; exit 1

# View products.
open: $(bldpdf); $(open) $^
openpngs: $(bldpngs); $(open) $^
openpdfs: $(bldpdfs) $(bldpdf); $(open) $^
openall: openpngs openpdfs

# Drop products.
clean cleantmp cleanone cleanpngs cleanpdfs cleandocs cleanall: clean%: ; $(rm) -rf $(bld$(or $*,tmp))

# Default target.
.DEFAULT_GOAL := build

# Non-file targets.
.PHONY: build t test check open openpngs openpdfs openall clean cleantmp cleanone cleanpngs cleanpdfs cleandocs cleanall
