set(SUBDIRS minpack) #odepack minpack colnew quadpack toms715 dfftpack slatec)
foreach (d ${SUBDIRS})
  add_subdirectory(${d})
endforeach (d)
