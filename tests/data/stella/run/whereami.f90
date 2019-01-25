PROGRAM whereami
  CHARACTER(len=255) :: cwd
  CALL getcwd(cwd)
  WRITE(*,'(/A, A)') 'I am in ', TRIM(cwd)
END PROGRAM whereami