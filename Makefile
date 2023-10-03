#!/usr/bin/make -f

# Default group.
groups := A15
g := $(or $g,$(firstword $(groups)))

# Default version.
versions := draft final
v := $(or $v,$(firstword $(versions)))

# Call make (again).
make = -r -j 4 -C $1 v=$2 $3
again = printf -v ident -- '%s %s %s' $1 $2 $3 \
				&& trap 'printf -- "\n$$ident: exit $$?.\n"' EXIT HUP INT TERM \
				&& printf -- "$$ident: exec \$$(MAKE) $(make)\n\n---\n" \
				&& TIMEFORMAT=$$'--- (%2lR elapsed %P%% busy)' \
				&& time $(MAKE) $(make)

# Update docs from here.
.DEFAULT_GOAL := docs
docs: $(groups);
	@set -eux && for g in $^; do cp -va $$g/build/$v/{$g.pdf,*.png} docs/$g/ && cp -va $$g/build/$v/$g.yml docs/_data/; done

# Redirect everything else.
$(groups): ; @$(call again,$@,$v,build)
$(versions): ; @$(call again,$g,$@,build)
.DEFAULT: ; @$(call again,$g,$v,$@)
.PHONY: docs $(groups) $(versions)
.SUFFIXES:
