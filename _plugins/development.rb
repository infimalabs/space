# Avoid redirect loop.
# https://github.com/jekyll/jekyll/issues/6459
module JekyllRedirectFrom
  class Generator < Jekyll::Generator
    alias_method :_redirectable_document?, :redirectable_document?
    def redirectable_document?(doc) false end
  end
end
