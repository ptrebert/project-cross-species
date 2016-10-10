# Comments on annotation config files

## encode_ctype_match.json

Matches between human and mouse cell types
based on either designated "pairs"
(labeled as analogs by ENCODE) or "close enough" argument.
Pairs that are wrong like "kidney - GM12878" represent
intentional mismatches used as negative controls
where needed/appropriate in the code.