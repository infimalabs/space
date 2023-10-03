#!/usr/bin/make -f

# Default group.
groups := A15
g := $(or $g,$(firstword $(groups)))

# Default version.
versions := draft final
v := $(or $v,$(firstword $(versions)))

# Build tooling.
install = install

# Call make (again).
make = $(MAKE) -r -j 4 -C $1 v=$2 $3
again = printf -v ident -- '%s %s %s' $1 $2 $3 \
				&& trap 'printf -- "\n$$ident: exit $$?.\n"' EXIT HUP INT TERM \
				&& printf -- "$$ident: exec \$$(MAKE) $(make)\n\n---\n" \
				&& TIMEFORMAT=$$'--- (%2lR elapsed %P%% busy)' \
				&& time $(make)

# Update docs from here.
.DEFAULT_GOAL := docs
docs: $(groups); @for g in $^; do (set -eux \
	&& $(install) -v -p -S -m 0644 $$g/_build/$v/*.png $$g \
	&& $(install) -v -p -S -m 0644 $$g/_build/$v/$$g.{md,pdf} . \
	); done

# Redirect everything else.
$(groups): ; @$(call again,$@,$v,build)
$(versions): ; @$(call again,$g,$@,build)
.DEFAULT: ; @$(call again,$g,$v,$@)
.PHONY: docs $(groups) $(versions)
.SUFFIXES:
