-- Drop \pdfvariables.
function Para(elem) if elem.content:find(pandoc.Str("suppressoptionalinfo")) then return {} else return elem end end

-- Drop \begin{opening}.
function Div(elem) if elem.classes:find("opening") then return {} else return elem end end
