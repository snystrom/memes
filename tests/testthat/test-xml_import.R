
tt_xml <- dremeR:::duplicate_file("tt_drememerge_dev/tomtom.xml")

importTomTomXML(tt_xml)

# testthat NULL tomtom result returns empty columns
importTomTomXML("inst/extdata/tomtom_ex/tomtom_bad.xml")
