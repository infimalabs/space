-- Fix url assets.
function _assets (url, t) return url:sub(1, 1) == "/" and url or "/A15/"..url end

-- Update <img src>.
function Image (img) img.src = _assets(img.src); return img end

-- Drop \pdfvariables.
function Para(elem) if elem.content:find(pandoc.Str("suppressoptionalinfo")) then return {} else return elem end end

-- Drop \begin{opening}.
function Div(elem) if elem.classes:find("opening") then return {} else return elem end end
