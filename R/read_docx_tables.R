# Read all tables of a .docx in the accepted view of tracked changes:
# text of <w:t> (inserted text included, <w:delText> excluded); rows marked as deleted are dropped.
read_docx_tables <- function(path) {
  td <- tempfile(); utils::unzip(path, "word/document.xml", exdir = td)
  doc <- xml2::read_xml(file.path(td, "word/document.xml")); ns <- xml2::xml_ns(doc)
  tbls <- xml2::xml_find_all(doc, "//w:body/w:tbl", ns)
  purrr::map(tbls, function(t) {
    rows <- xml2::xml_find_all(t, "./w:tr", ns)
    rows <- rows[!purrr::map_lgl(rows, ~ length(xml2::xml_find_all(.x, "./w:trPr/w:del", ns)) > 0)]
    purrr::map(rows, function(r) purrr::map_chr(xml2::xml_find_all(r, "./w:tc", ns),
      ~ paste(xml2::xml_text(xml2::xml_find_all(.x, ".//w:t", ns)), collapse = "")))
  })
}

# Body paragraphs (outside tables) in the accepted view of tracked changes
read_docx_paragraphs <- function(path) {
  td <- tempfile(); utils::unzip(path, "word/document.xml", exdir = td)
  doc <- xml2::read_xml(file.path(td, "word/document.xml")); ns <- xml2::xml_ns(doc)
  ps <- xml2::xml_find_all(doc, "//w:body/w:p", ns)
  txt <- purrr::map_chr(ps, ~ paste(xml2::xml_text(xml2::xml_find_all(.x, ".//w:t", ns)), collapse = ""))
  txt[nzchar(trimws(txt))]
}
