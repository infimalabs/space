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

# Build dirs.
bldall = _build
bldone = $(bldall)/$v
bldtmp = $(bldone)/tmp

# Build products.
bldpdfs = $(bldone)/$g.pdf
bldmark = $(bldone)/index.md
bldrtfm = $(bldone)/README.md
bldpngs = $(addprefix $(bldone)/,$(pngs:.txt=))

# Build tooling.
install = install -v
latexmk = latexmk
pandoc = pandoc
python = python3
mkdir = mkdir -v
open = open
cp = cp -v
rm = rm -v

# Version commands.
finalstampv = false

# Default commands.
figc = SAVEFIG=$@ $(python) $(addprefix ./,$^)
pre = \newcommand\version{$v} \
	\newcommand\group{$g} \
	\newcommand\stamp{$(or $($(v)stampv),true)} \
	\newcommand\setctime{\DTMsavetimestamp{ctime}{$(gittop)}} \
	\newcommand\setmtime{\DTMsavetimestamp{mtime}{$(gitmod)}}
latexc = $(cp) -a $(docs) $(bldone) \
				 && $(latexmk) -cd -outdir=tmp -jobname=$(@F).tmp -usepretex='$(pre)' -pdflua -quiet -halt-on-error -Werror $(bldone)/$(<F) \
				 && $(cp) -a $(bldtmp)/$(@F).tmp.pdf $@
docsc = $(rm) -f $@ \
				&& $(pandoc) $g.tex --output=$@ --lua-filter=$g.lua --resource-path=$(bldtmp) \
				--fail-if-warnings --standalone --citeproc --default-image-extension=png \
				-M toc-title:Table\ of\ Contents -M date:$(gittop) -M last_modified_at:$(gitmod) $1
htmlc = $(call docsc,--mathml --to html5+gfm_auto_identifiers --template=$g.tpl --toc --ascii --wrap=none)
markc = $(call docsc,--mathjax --to gfm-yaml_metadata_block)

# Ensure dirs.
$(bldone) $(bldtmp): ; $(mkdir) -p $@

# Compile figures.
$(bldone)/%.png: $g.py %.png.txt | $(bldone); $(figc)

# Compile content.
$(bldone)/%.pdf: $(docs) $(bldpngs) $(MAKEFILE_LIST) | $(bldtmp); $(latexc)

# Compile README.md.
$(bldrtfm): $(srcs) $(bldpngs) $(MAKEFILE_LIST); $(markc)

# Compile website.
$(bldmark): $(srcs) $(bldpngs) $(MAKEFILE_LIST) $(bldrtfm); $(htmlc)

# Group targets.
build: $(bldrtfm) $(bldmark) $(bldpdfs)

# Group tests.
t test check: build; @echo 'TODO($<): $@'

# Install products.
install: build test; $(install) -p -m 0644 $(bldpngs) $(bldrtfm) $(bldmark) $(bldpdfs) $(or $(DESTDIR),.)

# View products.
openpdfs open: $(bldpdfs); $(open) $^
openpngs: $(bldpngs); $(open) $^
openall: openpngs openpdfs

# Drop products.
clean: cleanpdfs cleantmp
cleantmp cleanone cleanpngs cleanpdfs cleanmark cleanrtfm cleanall: clean%: ; $(rm) -rf $(bld$*)

# Default target.
.DEFAULT_GOAL := install

# Non-file targets.
.PHONY: install build t test check open openpngs openpdfs openall \
	clean cleantmp cleanone cleanpngs cleanpdfs cleanmark cleanrtfm cleanall
