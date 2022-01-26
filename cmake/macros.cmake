# This macro expects a target (TARGT) and scope (SCOPE),
# based on which it will connect the sources to it.
# The list of pre-initialized source file names should be passed
# as an unnamed argument after the named arguments.
# Note that in macros, variables are in their string representations.
#    TARGT: TARGET (library, executable etc.)
#    SCOPE: PRIVATE | PUBLIC | INTERFACE
#    ARGN: List of unnamed keywords and arguments
macro(multi_src_target TARGT SCOPE)
	foreach(idx IN LISTS ${ARGN})
		target_sources(${TARGT} ${SCOPE} ${idx})
     	endforeach()
endmacro()
