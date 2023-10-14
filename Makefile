#!/usr/bin/make -f

# Default group.
groups := A15
g := $(or $g,$(firstword $(groups)))

# Default version.
versions := draft final
v := $(or $v,$(firstword $(versions)))

# Build tooling.
mv = mv -v

# Call make (again).
make = -r -j 4 -C $1 v=$2 $3
again = printf -v ident -- '%s %s %s' $1 $2 $3 \
				&& trap 'printf -- "\n$$ident: exit $$?.\n"' EXIT HUP INT TERM \
				&& printf -- "$$ident: exec \$$(MAKE) $(make)\n\n---\n" \
				&& TIMEFORMAT=$$'--- (%2lR elapsed %P%% busy)' \
				&& time $(MAKE) $(make)

# Install PDFs here.
.DEFAULT_GOAL := install
install: $(groups); @$(mv) $(groups:=/*.pdf) .

# Redirect everything else.
$(groups): ; @$(call again,$@,$v,install)
$(versions): ; @$(call again,$g,$@,install)
.DEFAULT: ; @$(call again,$g,$v,$@)
.PHONY: install $(groups) $(versions)
.SUFFIXES:
