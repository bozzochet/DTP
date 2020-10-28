
configure_file(weightfield2.root weightfield2.root COPYONLY)

install(
  FILES
    weightfield2.root
  DESTINATION
    "${INSTALL_DATA_DIR}"
)
