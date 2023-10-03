#!/usr/bin/make -f

# Group/Version.
override g := A15
override v := $(or $v,draft)

# Group sources.
docs := $g.tex $g.bib
pngs := $(notdir $(wildcard *.png.txt))
srcs := $g.tpl $g.lua $g.py $(pngs) $(docs)

# Git constants.
gittop := $(shell TZ=UTC0 git log -1 --date=format-local:%Y-%m-%dT%H:%M:%SZ --format=%ad --max-parents=0)
gitmod := $(shell TZ=UTC0 git log -1 --date=format-local:%Y-%m-%dT%H:%M:%SZ --format=%ad -- $(docs) $(pngs))
giturl := $(subst git@,https://,$(subst :,/,$(subst .git,,$(shell git config remote.origin.url))))

# Makefile utils.
reverse = $(if $(1),$(call reverse,$(wordlist 2,$(words $(1)),$(1)))) $(firstword $(1))

# Build dirs.
bldall = _build
bldone = $(bldall)/$v
bldtmp = $(bldone)/tmp

# Build products.
blddocs = $(bldone)/$g.md
bldpdfs = $(bldone)/$g.pdf
bldpngs = $(addprefix $(bldone)/,$(pngs:.txt=))

# Build tooling.
latexmk = latexmk
pandoc = pandoc
python = python3
mkdir = mkdir -v
touch = touch
open = open
cp = cp -v
mv = mv -v
rm = rm -v

# Version commands.
finalstampv = false

# Default commands.
figc = SAVEFIG=$@ $(python) $(addprefix ./,$^)
pre1 = \newcommand\version{$v}\newcommand\group{$g}\newcommand\stamp{$(or $($(v)stampv),true)}
pre2 = \newcommand\setgraphicspath{\graphicspath{$(addprefix {,$(addsuffix },$(call reverse,$(sort $(^D)))))}}
latexc = $(latexmk) -outdir=$(bldtmp) -jobname=$(@F).tmp -usepretex='$(pre1)$(pre2)' -pdflua -quiet -halt-on-error -Werror $< \
				 && $(cp) -a $(bldtmp)/$(@F).tmp.pdf $@.tmp.pdf \
				 && $(touch) -c -d $(gitmod) $@.tmp.pdf \
				 && $(mv) $@{.tmp.pdf,}
docsc = $(rm) -f $@ \
				&& $(pandoc) $g.tex --output=$@ --lua-filter=$g.lua --template=$g.tpl --resource-path=$(bldtmp) --fail-if-warnings --toc \
				--standalone --wrap=none --ascii --mathml --citeproc --default-image-extension=png --to=html5+gfm_auto_identifiers \
				-M toc-title:Table\ of\ Contents -M date:$(gittop) -M last_modified_at:$(gitmod) -M g:$g -M v:$v

# Ensure dirs.
$(bldone) $(bldtmp): ; $(mkdir) -p $@

# Compile figures.
$(bldone)/%.png: $g.py %.png.txt | $(bldone); $(figc)

# Compile content.
$(bldone)/%.pdf: %.tex $(bldpngs) $(MAKEFILE_LIST) | $(bldtmp); $(latexc)

# Compile website.
$(blddocs): $(srcs) $(bldpngs) $(MAKEFILE_LIST); $(docsc)

# Group targets.
build: $(blddocs) $(bldpdfs)

# Group tests.
t test check: build; @echo 'TODO($<): $@'; exit 1

# View products.
openpdfs open: $(bldpdfs); $(open) $^
openpngs: $(bldpngs); $(open) $^
openall: openpngs openpdfs

# Drop products.
clean: cleanpdfs cleantmp
cleantmp cleanone cleanpngs cleanpdfs cleandocs cleanall: clean%: ; $(rm) -rf $(bld$(or $*,tmp))

# Default target.
.DEFAULT_GOAL := build

# Non-file targets.
.PHONY: build t test check open openpngs openpdfs openall clean cleantmp cleanone cleanpngs cleanpdfs cleandocs cleanall
