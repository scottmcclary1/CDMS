!
!*** $Revision: 1.3.2.4 $
!*** $Date: 2007/03/22 20:40:24 $
!***
!*** Copyright 1985-2007 Intel Corporation.  All Rights Reserved.
!***
!*** The source code contained or described herein and all documents related
!*** to the source code ("Material") are owned by Intel Corporation or its
!*** suppliers or licensors.  Title to the Material remains with Intel
!*** Corporation or its suppliers and licensors.  The Material is protected
!*** by worldwide copyright laws and treaty provisions.  No part of the
!*** Material may be used, copied, reproduced, modified, published, uploaded,
!*** posted, transmitted, distributed, or disclosed in any way without
!*** Intel's prior express written permission.
!***
!*** No license under any patent, copyright, trade secret or other
!*** intellectual property right is granted to or conferred upon you by
!*** disclosure or delivery of the Materials, either expressly, by
!*** implication, inducement, estoppel or otherwise.  Any license under such
!*** intellectual property rights must be express and approved by Intel in
!*** writing.
!***
!*** Portions of this software are protected under the following patents:
!***     U.S. Patent 5,812,852
!***     U.S. Patent 6,792,599
!***

!***
!*** Some of the directives for the following routine extend past column 72,
!*** so process this file in 132-column mode.
!***

!dec$ fixedformlinesize:132

        module omp_lib_kinds
        integer omp_integer_kind
        integer omp_logical_kind
        integer omp_real_kind
        integer omp_lock_kind
        integer omp_nest_lock_kind
        integer openmp_version
        integer kmp_pointer_kind
        integer kmp_size_t_kind
        parameter (omp_integer_kind = 4)
        parameter (omp_logical_kind = 4)
        parameter (omp_real_kind = 4)
        parameter (omp_lock_kind = int_ptr_kind())
        parameter (omp_nest_lock_kind = int_ptr_kind())
        parameter (openmp_version = 200505)
        parameter (kmp_pointer_kind = int_ptr_kind())
        parameter (kmp_size_t_kind = int_ptr_kind())
        end module omp_lib_kinds
        module omp_lib
          use omp_lib_kinds

!***
!*** omp_* entry points
!***

          interface
            subroutine omp_set_num_threads(nthreads)
              use omp_lib_kinds
              integer (kind=omp_integer_kind) nthreads
            end subroutine omp_set_num_threads
          end interface

          interface
            subroutine omp_set_dynamic(enable)
              use omp_lib_kinds
              logical (kind=omp_logical_kind) enable
            end subroutine omp_set_dynamic
          end interface

          interface
            subroutine omp_set_nested(enable)
              use omp_lib_kinds
              logical (kind=omp_logical_kind) enable
            end subroutine omp_set_nested
          end interface

          interface
            function omp_get_num_threads()
              use omp_lib_kinds
              integer (kind=omp_integer_kind) omp_get_num_threads
            end function omp_get_num_threads
          end interface

          interface
            function omp_get_max_threads()
              use omp_lib_kinds
              integer (kind=omp_integer_kind) omp_get_max_threads
            end function omp_get_max_threads
          end interface

          interface
            function omp_get_thread_num()
              use omp_lib_kinds
              integer (kind=omp_integer_kind) omp_get_thread_num
            end function omp_get_thread_num
          end interface

          interface
            function omp_get_num_procs()
              use omp_lib_kinds
              integer (kind=omp_integer_kind) omp_get_num_procs
            end function omp_get_num_procs
          end interface

          interface
            function omp_in_parallel()
              use omp_lib_kinds
              logical (kind=omp_logical_kind) omp_in_parallel
            end function omp_in_parallel
          end interface

          interface
            function omp_get_dynamic()
              use omp_lib_kinds
              logical (kind=omp_logical_kind) omp_get_dynamic
            end function omp_get_dynamic
          end interface

          interface
            function omp_get_nested()
              use omp_lib_kinds
              logical (kind=omp_logical_kind) omp_get_nested
            end function omp_get_nested
          end interface

          interface
            function omp_get_wtime()
              use omp_lib_kinds
              double precision omp_get_wtime
            end function omp_get_wtime
          end interface

          interface
            function omp_get_wtick ()
              use omp_lib_kinds
              double precision omp_get_wtick
            end function omp_get_wtick
          end interface

          interface
            subroutine omp_init_lock(lockvar)
              use omp_lib_kinds
              integer (kind=omp_lock_kind) lockvar
            end subroutine omp_init_lock
          end interface

          interface
            subroutine omp_destroy_lock(lockvar)
              use omp_lib_kinds
              integer (kind=omp_lock_kind) lockvar
            end subroutine omp_destroy_lock
          end interface

          interface
            subroutine omp_set_lock(lockvar)
              use omp_lib_kinds
              integer (kind=omp_lock_kind) lockvar
            end subroutine omp_set_lock
          end interface

          interface
            subroutine omp_unset_lock(lockvar)
              use omp_lib_kinds
              integer (kind=omp_lock_kind) lockvar
            end subroutine omp_unset_lock
          end interface

          interface
            function omp_test_lock(lockvar)
              use omp_lib_kinds
              logical (kind=omp_logical_kind) omp_test_lock
              integer (kind=omp_lock_kind) lockvar
            end function omp_test_lock
          end interface

          interface
            subroutine omp_init_nest_lock(lockvar)
              use omp_lib_kinds
              integer (kind=omp_nest_lock_kind) lockvar
            end subroutine omp_init_nest_lock
          end interface

          interface
            subroutine omp_destroy_nest_lock(lockvar)
              use omp_lib_kinds
              integer (kind=omp_nest_lock_kind) lockvar
            end subroutine omp_destroy_nest_lock
          end interface

          interface
            subroutine omp_set_nest_lock(lockvar)
              use omp_lib_kinds
              integer (kind=omp_nest_lock_kind) lockvar
            end subroutine omp_set_nest_lock
          end interface

          interface
            subroutine omp_unset_nest_lock(lockvar)
              use omp_lib_kinds
              integer (kind=omp_nest_lock_kind) lockvar
            end subroutine omp_unset_nest_lock
          end interface

          interface
            function omp_test_nest_lock(lockvar)
              use omp_lib_kinds
              integer (kind=omp_integer_kind) omp_test_nest_lock
              integer (kind=omp_nest_lock_kind) lockvar
            end function omp_test_nest_lock
          end interface

!***
!*** kmp_* entry points
!***

          interface
            subroutine kmp_set_parallel_name(name)
              use omp_lib_kinds
              character*(*) name
            end subroutine kmp_set_parallel_name
          end interface

          interface
            subroutine kmp_set_stacksize(size)
              use omp_lib_kinds
              integer (kind=omp_integer_kind) size
            end subroutine kmp_set_stacksize
          end interface

          interface
            subroutine kmp_set_stacksize_s(size)
              use omp_lib_kinds
              integer (kind=kmp_size_t_kind) size
            end subroutine kmp_set_stacksize_s
          end interface

          interface
            subroutine kmp_set_blocktime(msec)
              use omp_lib_kinds
              integer (kind=omp_integer_kind) msec
            end subroutine kmp_set_blocktime
          end interface

          interface
            subroutine kmp_set_library_serial()
              use omp_lib_kinds
            end subroutine kmp_set_library_serial
          end interface

          interface
            subroutine kmp_set_library_turnaround()
              use omp_lib_kinds
            end subroutine kmp_set_library_turnaround
          end interface

          interface
            subroutine kmp_set_library_throughput()
              use omp_lib_kinds
            end subroutine kmp_set_library_throughput
          end interface

          interface
            subroutine kmp_set_library(libnum)
              use omp_lib_kinds
              integer (kind=omp_integer_kind) libnum
            end subroutine kmp_set_library
          end interface

          interface
            subroutine kmp_set_stats(enable)
              use omp_lib_kinds
              logical (kind=omp_logical_kind) enable
            end subroutine kmp_set_stats
          end interface

          interface
            function kmp_get_stacksize()
              use omp_lib_kinds
              integer (kind=omp_integer_kind) kmp_get_stacksize
            end function kmp_get_stacksize
          end interface

          interface
            function kmp_get_stacksize_s()
              use omp_lib_kinds
              integer (kind=kmp_size_t_kind) kmp_get_stacksize_s
            end function kmp_get_stacksize_s
          end interface

          interface
            function kmp_get_blocktime()
              use omp_lib_kinds
              integer (kind=omp_integer_kind) kmp_get_blocktime
            end function kmp_get_blocktime
          end interface

          interface
            function kmp_get_library()
              use omp_lib_kinds
              integer (kind=omp_integer_kind) kmp_get_library
            end function kmp_get_library
          end interface

          interface
            function kmp_malloc(size)
              use omp_lib_kinds
              integer (kind=kmp_pointer_kind) kmp_malloc
              integer (kind=kmp_size_t_kind) size
            end function kmp_malloc
          end interface

          interface
            function kmp_calloc(nelem, elsize)
              use omp_lib_kinds
              integer (kind=kmp_pointer_kind) kmp_calloc
              integer (kind=kmp_size_t_kind) nelem
              integer (kind=kmp_size_t_kind) elsize
            end function kmp_calloc
          end interface

          interface
            function kmp_realloc(ptr, size)
              use omp_lib_kinds
              integer (kind=kmp_pointer_kind) kmp_realloc
              integer (kind=kmp_pointer_kind) ptr
              integer (kind=kmp_size_t_kind) size
            end function kmp_realloc
          end interface

          interface
            subroutine kmp_free(ptr)
              use omp_lib_kinds
              integer (kind=kmp_pointer_kind) ptr
            end subroutine kmp_free
          end interface

!***
!*** kmp_* entry points for Cluster OpenMP
!***

          interface
            function kmp_sharable_malloc(size)
              use omp_lib_kinds
              integer (kind=kmp_pointer_kind) kmp_sharable_malloc
              integer (kind=kmp_size_t_kind) size
            end function kmp_sharable_malloc
          end interface

          interface
            function kmp_aligned_sharable_malloc(size)
              use omp_lib_kinds
              integer (kind=kmp_pointer_kind) kmp_aligned_sharable_malloc
              integer (kind=kmp_size_t_kind) size
            end function kmp_aligned_sharable_malloc
          end interface

          interface
            function kmp_sharable_calloc(nelem, elsize)
              use omp_lib_kinds
              integer (kind=kmp_pointer_kind) kmp_sharable_calloc
              integer (kind=kmp_size_t_kind) nelem
              integer (kind=kmp_size_t_kind) elsize
            end function kmp_sharable_calloc
          end interface

          interface
            function kmp_sharable_realloc(ptr, size)
              use omp_lib_kinds
              integer (kind=kmp_pointer_kind) kmp_sharable_realloc
              integer (kind=kmp_pointer_kind) ptr
              integer (kind=kmp_size_t_kind) size
            end function kmp_sharable_realloc
          end interface

          interface
            subroutine kmp_sharable_free(ptr)
              use omp_lib_kinds
              integer (kind=kmp_pointer_kind) ptr
            end subroutine kmp_sharable_free
          end interface

          interface
            subroutine kmp_lock_cond_wait(lockvar)
              use omp_lib_kinds
              integer (kind=omp_lock_kind) lockvar
            end subroutine kmp_lock_cond_wait
          end interface

          interface
            subroutine kmp_lock_cond_signal(lockvar)
              use omp_lib_kinds
              integer (kind=omp_lock_kind) lockvar
            end subroutine kmp_lock_cond_signal
          end interface

          interface
            subroutine kmp_lock_cond_broadcast(lockvar)
              use omp_lib_kinds
              integer (kind=omp_lock_kind) lockvar
            end subroutine kmp_lock_cond_broadcast
          end interface

          interface
            subroutine kmp_nest_lock_cond_wait(lockvar)
              use omp_lib_kinds
              integer (kind=omp_nest_lock_kind) lockvar
            end subroutine kmp_nest_lock_cond_wait
          end interface

          interface
            subroutine kmp_nest_lock_cond_signal(lockvar)
              use omp_lib_kinds
              integer (kind=omp_nest_lock_kind) lockvar
            end subroutine kmp_nest_lock_cond_signal
          end interface

          interface
            subroutine kmp_nest_lock_cond_broadcast(lockvar)
              use omp_lib_kinds
              integer (kind=omp_nest_lock_kind) lockvar
            end subroutine kmp_nest_lock_cond_broadcast
          end interface

          interface
            function kmp_get_num_processes()
              use omp_lib_kinds
              integer (kind=omp_integer_kind) kmp_get_num_processes
            end function kmp_get_num_processes
          end interface

          interface
            function kmp_get_process_num()
              use omp_lib_kinds
              integer (kind=omp_integer_kind) kmp_get_process_num
            end function kmp_get_process_num
          end interface

          interface
            function kmp_get_process_thread_num()
              use omp_lib_kinds
              integer (kind=omp_integer_kind) kmp_get_process_thread_num
            end function kmp_get_process_thread_num
          end interface

          interface
            subroutine kmp_set_warnings_on()
              use omp_lib_kinds
            end subroutine kmp_set_warnings_on
          end interface

          interface
            subroutine kmp_set_warnings_off()
              use omp_lib_kinds
            end subroutine kmp_set_warnings_off
          end interface

          interface
            function kmp_is_sharable(ptr)
              use omp_lib_kinds
              logical (kind=omp_logical_kind) kmp_is_sharable
              integer (kind=kmp_pointer_kind) ptr
            end function kmp_is_sharable
          end interface

          interface
            subroutine kmp_deferred_atomic_add_i4(addr, val)
              integer(kind=4) addr
              integer(kind=4) val
            end subroutine kmp_deferred_atomic_add_i4
          end interface

          interface
            subroutine kmp_deferred_atomic_add_i8(addr, val)
              integer(kind=8) addr
              integer(kind=8) val
            end subroutine kmp_deferred_atomic_add_i8
          end interface

          interface
            subroutine kmp_deferred_atomic_add_r4(addr, val)
              real(kind=4) addr
              real(kind=4) val
            end subroutine kmp_deferred_atomic_add_r4
          end interface

          interface
            subroutine kmp_deferred_atomic_add_r8(addr, val)
              real(kind=8) addr
              real(kind=8) val
            end subroutine kmp_deferred_atomic_add_r8
          end interface

!dec$ if defined(_WIN32)
!dec$   if defined(_WIN64) .or. defined(_M_IA64) .or. defined(_M_AMD64)

!***
!*** The Fortran entry points must be in uppercase, even if the /Qlowercase
!*** option is specified.  The alias attribute ensures that the specified
!*** string is used as the entry point.
!***
!*** On the Windows IA-32 architecture, the Fortran entry points have an
!*** underscore prepended.  On the Windows Intel(R) 64 and Intel(R) Itanium(R)
!*** architectures, no underscore is prepended.
!***

!dec$ attributes alias:'OMP_SET_NUM_THREADS' :: omp_set_num_threads
!dec$ attributes alias:'OMP_SET_DYNAMIC' :: omp_set_dynamic
!dec$ attributes alias:'OMP_SET_NESTED' :: omp_set_nested
!dec$ attributes alias:'OMP_GET_NUM_THREADS' :: omp_get_num_threads
!dec$ attributes alias:'OMP_GET_MAX_THREADS' :: omp_get_max_threads
!dec$ attributes alias:'OMP_GET_THREAD_NUM' :: omp_get_thread_num
!dec$ attributes alias:'OMP_GET_NUM_PROCS' :: omp_get_num_procs
!dec$ attributes alias:'OMP_IN_PARALLEL' :: omp_in_parallel
!dec$ attributes alias:'OMP_GET_DYNAMIC' :: omp_get_dynamic
!dec$ attributes alias:'OMP_GET_NESTED' :: omp_get_nested
!dec$ attributes alias:'OMP_GET_WTIME' :: omp_get_wtime
!dec$ attributes alias:'OMP_GET_WTICK' :: omp_get_wtick

!dec$ attributes alias:'omp_init_lock' :: omp_init_lock
!dec$ attributes alias:'omp_destroy_lock' :: omp_destroy_lock
!dec$ attributes alias:'omp_set_lock' :: omp_set_lock
!dec$ attributes alias:'omp_unset_lock' :: omp_unset_lock
!dec$ attributes alias:'omp_test_lock' :: omp_test_lock
!dec$ attributes alias:'omp_init_nest_lock' :: omp_init_nest_lock
!dec$ attributes alias:'omp_destroy_nest_lock' :: omp_destroy_nest_lock
!dec$ attributes alias:'omp_set_nest_lock' :: omp_set_nest_lock
!dec$ attributes alias:'omp_unset_nest_lock' :: omp_unset_nest_lock
!dec$ attributes alias:'omp_test_nest_lock' :: omp_test_nest_lock

!dec$ attributes alias:'KMP_SET_PARALLEL_NAME'::kmp_set_parallel_name
!dec$ attributes alias:'KMP_SET_STACKSIZE'::kmp_set_stacksize
!dec$ attributes alias:'KMP_SET_STACKSIZE_S'::kmp_set_stacksize_s
!dec$ attributes alias:'KMP_SET_BLOCKTIME'::kmp_set_blocktime
!dec$ attributes alias:'KMP_SET_LIBRARY_SERIAL'::kmp_set_library_serial
!dec$ attributes alias:'KMP_SET_LIBRARY_TURNAROUND'::kmp_set_library_turnaround
!dec$ attributes alias:'KMP_SET_LIBRARY_THROUGHPUT'::kmp_set_library_throughput
!dec$ attributes alias:'KMP_SET_LIBRARY'::kmp_set_library
!dec$ attributes alias:'KMP_SET_STATS'::kmp_set_stats
!dec$ attributes alias:'KMP_GET_STACKSIZE'::kmp_get_stacksize
!dec$ attributes alias:'KMP_GET_STACKSIZE_S'::kmp_get_stacksize_s
!dec$ attributes alias:'KMP_GET_BLOCKTIME'::kmp_get_blocktime
!dec$ attributes alias:'KMP_GET_LIBRARY'::kmp_get_library
!dec$ attributes alias:'KMP_MALLOC'::kmp_malloc
!dec$ attributes alias:'KMP_CALLOC'::kmp_calloc
!dec$ attributes alias:'KMP_REALLOC'::kmp_realloc
!dec$ attributes alias:'KMP_FREE'::kmp_free

!dec$ attributes alias:'KMP_SHARABLE_MALLOC'::kmp_sharable_malloc
!dec$ attributes alias:'KMP_ALIGNED_SHARABLE_MALLOC'::kmp_aligned_sharable_malloc
!dec$ attributes alias:'KMP_SHARABLE_CALLOC'::kmp_sharable_calloc
!dec$ attributes alias:'KMP_SHARABLE_REALLOC'::kmp_sharable_realloc
!dec$ attributes alias:'KMP_SHARABLE_FREE'::kmp_sharable_free
!dec$ attributes alias:'KMP_LOCK_COND_WAIT'::kmp_lock_cond_wait
!dec$ attributes alias:'KMP_LOCK_COND_SIGNAL'::kmp_lock_cond_signal
!dec$ attributes alias:'KMP_LOCK_COND_BROADCAST'::kmp_lock_cond_broadcast
!dec$ attributes alias:'KMP_NEST_LOCK_COND_WAIT'::kmp_nest_lock_cond_wait
!dec$ attributes alias:'KMP_NEST_LOCK_COND_SIGNAL'::kmp_nest_lock_cond_signal
!dec$ attributes alias:'KMP_NEST_LOCK_COND_BROADCAST'::kmp_nest_lock_cond_broadcast
!dec$ attributes alias:'KMP_GET_NUM_PROCESSES'::kmp_get_num_processes
!dec$ attributes alias:'KMP_GET_PROCESS_NUM'::kmp_get_process_num
!dec$ attributes alias:'KMP_GET_PROCESS_THREAD_NUM'::kmp_get_process_thread_num
!dec$ attributes alias:'KMP_SET_WARNINGS_ON'::kmp_set_warnings_on
!dec$ attributes alias:'KMP_SET_WARNINGS_OFF'::kmp_set_warnings_off
!dec$ attributes alias:'KMP_IS_SHARABLE'::kmp_is_sharable
!dec$ attributes alias:'KMP_DEFERRED_ATOMIC_ADD_I4'::kmp_deferred_atomic_add_i4
!dec$ attributes alias:'KMP_DEFERRED_ATOMIC_ADD_I8'::kmp_deferred_atomic_add_i8
!dec$ attributes alias:'KMP_DEFERRED_ATOMIC_ADD_R4'::kmp_deferred_atomic_add_r4
!dec$ attributes alias:'KMP_DEFERRED_ATOMIC_ADD_R8'::kmp_deferred_atomic_add_r8

!dec$   else

!***
!*** On Windows IA32, the Fortran entry points have an underscore prepended.
!***

!dec$ attributes alias:'_OMP_SET_NUM_THREADS' :: omp_set_num_threads
!dec$ attributes alias:'_OMP_SET_DYNAMIC' :: omp_set_dynamic
!dec$ attributes alias:'_OMP_SET_NESTED' :: omp_set_nested
!dec$ attributes alias:'_OMP_GET_NUM_THREADS' :: omp_get_num_threads
!dec$ attributes alias:'_OMP_GET_MAX_THREADS' :: omp_get_max_threads
!dec$ attributes alias:'_OMP_GET_THREAD_NUM' :: omp_get_thread_num
!dec$ attributes alias:'_OMP_GET_NUM_PROCS' :: omp_get_num_procs
!dec$ attributes alias:'_OMP_IN_PARALLEL' :: omp_in_parallel
!dec$ attributes alias:'_OMP_GET_DYNAMIC' :: omp_get_dynamic
!dec$ attributes alias:'_OMP_GET_NESTED' :: omp_get_nested
!dec$ attributes alias:'_OMP_GET_WTIME' :: omp_get_wtime
!dec$ attributes alias:'_OMP_GET_WTICK' :: omp_get_wtick

!dec$ attributes alias:'_omp_init_lock' :: omp_init_lock
!dec$ attributes alias:'_omp_destroy_lock' :: omp_destroy_lock
!dec$ attributes alias:'_omp_set_lock' :: omp_set_lock
!dec$ attributes alias:'_omp_unset_lock' :: omp_unset_lock
!dec$ attributes alias:'_omp_test_lock' :: omp_test_lock
!dec$ attributes alias:'_omp_init_nest_lock' :: omp_init_nest_lock
!dec$ attributes alias:'_omp_destroy_nest_lock' :: omp_destroy_nest_lock
!dec$ attributes alias:'_omp_set_nest_lock' :: omp_set_nest_lock
!dec$ attributes alias:'_omp_unset_nest_lock' :: omp_unset_nest_lock
!dec$ attributes alias:'_omp_test_nest_lock' :: omp_test_nest_lock

!dec$ attributes alias:'_KMP_SET_PARALLEL_NAME'::kmp_set_parallel_name
!dec$ attributes alias:'_KMP_SET_STACKSIZE'::kmp_set_stacksize
!dec$ attributes alias:'_KMP_SET_STACKSIZE_S'::kmp_set_stacksize_s
!dec$ attributes alias:'_KMP_SET_BLOCKTIME'::kmp_set_blocktime
!dec$ attributes alias:'_KMP_SET_LIBRARY_SERIAL'::kmp_set_library_serial
!dec$ attributes alias:'_KMP_SET_LIBRARY_TURNAROUND'::kmp_set_library_turnaround
!dec$ attributes alias:'_KMP_SET_LIBRARY_THROUGHPUT'::kmp_set_library_throughput
!dec$ attributes alias:'_KMP_SET_LIBRARY'::kmp_set_library
!dec$ attributes alias:'_KMP_SET_STATS'::kmp_set_stats
!dec$ attributes alias:'_KMP_GET_STACKSIZE'::kmp_get_stacksize
!dec$ attributes alias:'_KMP_GET_STACKSIZE_S'::kmp_get_stacksize_s
!dec$ attributes alias:'_KMP_GET_BLOCKTIME'::kmp_get_blocktime
!dec$ attributes alias:'_KMP_GET_LIBRARY'::kmp_get_library
!dec$ attributes alias:'_KMP_MALLOC'::kmp_malloc
!dec$ attributes alias:'_KMP_CALLOC'::kmp_calloc
!dec$ attributes alias:'_KMP_REALLOC'::kmp_realloc
!dec$ attributes alias:'_KMP_FREE'::kmp_free

!dec$ attributes alias:'_KMP_SHARABLE_MALLOC'::kmp_sharable_malloc
!dec$ attributes alias:'_KMP_ALIGNED_SHARABLE_MALLOC'::kmp_aligned_sharable_malloc
!dec$ attributes alias:'_KMP_SHARABLE_CALLOC'::kmp_sharable_calloc
!dec$ attributes alias:'_KMP_SHARABLE_REALLOC'::kmp_sharable_realloc
!dec$ attributes alias:'_KMP_SHARABLE_FREE'::kmp_sharable_free
!dec$ attributes alias:'_KMP_LOCK_COND_WAIT'::kmp_lock_cond_wait
!dec$ attributes alias:'_KMP_LOCK_COND_SIGNAL'::kmp_lock_cond_signal
!dec$ attributes alias:'_KMP_LOCK_COND_BROADCAST'::kmp_lock_cond_broadcast
!dec$ attributes alias:'_KMP_NEST_LOCK_COND_WAIT'::kmp_nest_lock_cond_wait
!dec$ attributes alias:'_KMP_NEST_LOCK_COND_SIGNAL'::kmp_nest_lock_cond_signal
!dec$ attributes alias:'_KMP_NEST_LOCK_COND_BROADCAST'::kmp_nest_lock_cond_broadcast
!dec$ attributes alias:'_KMP_GET_NUM_PROCESSES'::kmp_get_num_processes
!dec$ attributes alias:'_KMP_GET_PROCESS_NUM'::kmp_get_process_num
!dec$ attributes alias:'_KMP_GET_PROCESS_THREAD_NUM'::kmp_get_process_thread_num
!dec$ attributes alias:'_KMP_SET_WARNINGS_ON'::kmp_set_warnings_on
!dec$ attributes alias:'_KMP_SET_WARNINGS_OFF'::kmp_set_warnings_off
!dec$ attributes alias:'_KMP_IS_SHARABLE'::kmp_is_sharable
!dec$ attributes alias:'_KMP_DEFERRED_ATOMIC_ADD_I4'::kmp_deferred_atomic_add_i4
!dec$ attributes alias:'_KMP_DEFERRED_ATOMIC_ADD_I8'::kmp_deferred_atomic_add_i8
!dec$ attributes alias:'_KMP_DEFERRED_ATOMIC_ADD_R4'::kmp_deferred_atomic_add_r4
!dec$ attributes alias:'_KMP_DEFERRED_ATOMIC_ADD_R8'::kmp_deferred_atomic_add_r8

!dec$   endif
!dec$ endif

!dec$ if defined(__linux)

!***
!*** The Linux entry points are in lowercase, with an underscore appended.
!***

!dec$ attributes alias:'omp_set_num_threads_'::omp_set_num_threads
!dec$ attributes alias:'omp_set_dynamic_'::omp_set_dynamic
!dec$ attributes alias:'omp_set_nested_'::omp_set_nested
!dec$ attributes alias:'omp_get_num_threads_'::omp_get_num_threads
!dec$ attributes alias:'omp_get_max_threads_'::omp_get_max_threads
!dec$ attributes alias:'omp_get_thread_num_'::omp_get_thread_num
!dec$ attributes alias:'omp_get_num_procs_'::omp_get_num_procs
!dec$ attributes alias:'omp_in_parallel_'::omp_in_parallel
!dec$ attributes alias:'omp_get_dynamic_'::omp_get_dynamic
!dec$ attributes alias:'omp_get_nested_'::omp_get_nested
!dec$ attributes alias:'omp_get_wtime_'::omp_get_wtime
!dec$ attributes alias:'omp_get_wtick_'::omp_get_wtick

!dec$ attributes alias:'omp_init_lock_'::omp_init_lock
!dec$ attributes alias:'omp_destroy_lock_'::omp_destroy_lock
!dec$ attributes alias:'omp_set_lock_'::omp_set_lock
!dec$ attributes alias:'omp_unset_lock_'::omp_unset_lock
!dec$ attributes alias:'omp_test_lock_'::omp_test_lock
!dec$ attributes alias:'omp_init_nest_lock_'::omp_init_nest_lock
!dec$ attributes alias:'omp_destroy_nest_lock_'::omp_destroy_nest_lock
!dec$ attributes alias:'omp_set_nest_lock_'::omp_set_nest_lock
!dec$ attributes alias:'omp_unset_nest_lock_'::omp_unset_nest_lock
!dec$ attributes alias:'omp_test_nest_lock_'::omp_test_nest_lock

!dec$ attributes alias:'kmp_set_parallel_name_'::kmp_set_parallel_name
!dec$ attributes alias:'kmp_set_stacksize_'::kmp_set_stacksize
!dec$ attributes alias:'kmp_set_stacksize_s_'::kmp_set_stacksize_s
!dec$ attributes alias:'kmp_set_blocktime_'::kmp_set_blocktime
!dec$ attributes alias:'kmp_set_library_serial_'::kmp_set_library_serial
!dec$ attributes alias:'kmp_set_library_turnaround_'::kmp_set_library_turnaround
!dec$ attributes alias:'kmp_set_library_throughput_'::kmp_set_library_throughput
!dec$ attributes alias:'kmp_set_library_'::kmp_set_library
!dec$ attributes alias:'kmp_set_stats_'::kmp_set_stats
!dec$ attributes alias:'kmp_get_stacksize_'::kmp_get_stacksize
!dec$ attributes alias:'kmp_get_stacksize_s_'::kmp_get_stacksize_s
!dec$ attributes alias:'kmp_get_blocktime_'::kmp_get_blocktime
!dec$ attributes alias:'kmp_get_library_'::kmp_get_library
!dec$ attributes alias:'kmp_malloc_'::kmp_malloc
!dec$ attributes alias:'kmp_calloc_'::kmp_calloc
!dec$ attributes alias:'kmp_realloc_'::kmp_realloc
!dec$ attributes alias:'kmp_free_'::kmp_free

!dec$ attributes alias:'kmp_sharable_malloc_'::kmp_sharable_malloc
!dec$ attributes alias:'kmp_aligned_sharable_malloc_'::kmp_aligned_sharable_malloc
!dec$ attributes alias:'kmp_sharable_calloc_'::kmp_sharable_calloc
!dec$ attributes alias:'kmp_sharable_realloc_'::kmp_sharable_realloc
!dec$ attributes alias:'kmp_sharable_free_'::kmp_sharable_free
!dec$ attributes alias:'kmp_lock_cond_wait_'::kmp_lock_cond_wait
!dec$ attributes alias:'kmp_lock_cond_signal_'::kmp_lock_cond_signal
!dec$ attributes alias:'kmp_lock_cond_broadcast_'::kmp_lock_cond_broadcast
!dec$ attributes alias:'kmp_nest_lock_cond_wait_'::kmp_nest_lock_cond_wait
!dec$ attributes alias:'kmp_nest_lock_cond_signal_'::kmp_nest_lock_cond_signal
!dec$ attributes alias:'kmp_nest_lock_cond_broadcast_'::kmp_nest_lock_cond_broadcast
!dec$ attributes alias:'kmp_get_num_processes_'::kmp_get_num_processes
!dec$ attributes alias:'kmp_get_process_num_'::kmp_get_process_num
!dec$ attributes alias:'kmp_get_process_thread_num_'::kmp_get_process_thread_num
!dec$ attributes alias:'kmp_set_warnings_on_'::kmp_set_warnings_on
!dec$ attributes alias:'kmp_set_warnings_off_'::kmp_set_warnings_off
!dec$ attributes alias:'kmp_is_sharable_'::kmp_is_sharable
!dec$ attributes alias:'kmp_deferred_atomic_add_i4_'::kmp_deferred_atomic_add_i4
!dec$ attributes alias:'kmp_deferred_atomic_add_i8_'::kmp_deferred_atomic_add_i8
!dec$ attributes alias:'kmp_deferred_atomic_add_r4_'::kmp_deferred_atomic_add_r4
!dec$ attributes alias:'kmp_deferred_atomic_add_r8_'::kmp_deferred_atomic_add_r8

!dec$ endif

!dec$ if defined(__APPLE__)

!***
!*** The Mac entry points are in lowercase, with an both an underscore
!*** appended and an underscore prepended.
!***

!dec$ attributes alias:'_omp_set_num_threads_'::omp_set_num_threads
!dec$ attributes alias:'_omp_set_dynamic_'::omp_set_dynamic
!dec$ attributes alias:'_omp_set_nested_'::omp_set_nested
!dec$ attributes alias:'_omp_get_num_threads_'::omp_get_num_threads
!dec$ attributes alias:'_omp_get_max_threads_'::omp_get_max_threads
!dec$ attributes alias:'_omp_get_thread_num_'::omp_get_thread_num
!dec$ attributes alias:'_omp_get_num_procs_'::omp_get_num_procs
!dec$ attributes alias:'_omp_in_parallel_'::omp_in_parallel
!dec$ attributes alias:'_omp_get_dynamic_'::omp_get_dynamic
!dec$ attributes alias:'_omp_get_nested_'::omp_get_nested
!dec$ attributes alias:'_omp_get_wtime_'::omp_get_wtime
!dec$ attributes alias:'_omp_get_wtick_'::omp_get_wtick

!dec$ attributes alias:'_omp_init_lock_'::omp_init_lock
!dec$ attributes alias:'_omp_destroy_lock_'::omp_destroy_lock
!dec$ attributes alias:'_omp_set_lock_'::omp_set_lock
!dec$ attributes alias:'_omp_unset_lock_'::omp_unset_lock
!dec$ attributes alias:'_omp_test_lock_'::omp_test_lock
!dec$ attributes alias:'_omp_init_nest_lock_'::omp_init_nest_lock
!dec$ attributes alias:'_omp_destroy_nest_lock_'::omp_destroy_nest_lock
!dec$ attributes alias:'_omp_set_nest_lock_'::omp_set_nest_lock
!dec$ attributes alias:'_omp_unset_nest_lock_'::omp_unset_nest_lock
!dec$ attributes alias:'_omp_test_nest_lock_'::omp_test_nest_lock

!dec$ attributes alias:'_kmp_set_parallel_name_'::kmp_set_parallel_name
!dec$ attributes alias:'_kmp_set_stacksize_'::kmp_set_stacksize
!dec$ attributes alias:'_kmp_set_stacksize_s_'::kmp_set_stacksize_s
!dec$ attributes alias:'_kmp_set_blocktime_'::kmp_set_blocktime
!dec$ attributes alias:'_kmp_set_library_serial_'::kmp_set_library_serial
!dec$ attributes alias:'_kmp_set_library_turnaround_'::kmp_set_library_turnaround
!dec$ attributes alias:'_kmp_set_library_throughput_'::kmp_set_library_throughput
!dec$ attributes alias:'_kmp_set_library_'::kmp_set_library
!dec$ attributes alias:'_kmp_set_stats_'::kmp_set_stats
!dec$ attributes alias:'_kmp_get_stacksize_'::kmp_get_stacksize
!dec$ attributes alias:'_kmp_get_stacksize_s_'::kmp_get_stacksize_s
!dec$ attributes alias:'_kmp_get_blocktime_'::kmp_get_blocktime
!dec$ attributes alias:'_kmp_get_library_'::kmp_get_library
!dec$ attributes alias:'_kmp_malloc_'::kmp_malloc
!dec$ attributes alias:'_kmp_calloc_'::kmp_calloc
!dec$ attributes alias:'_kmp_realloc_'::kmp_realloc
!dec$ attributes alias:'_kmp_free_'::kmp_free

!dec$ attributes alias:'_kmp_sharable_malloc_'::kmp_sharable_malloc
!dec$ attributes alias:'_kmp_aligned_sharable_malloc_'::kmp_aligned_sharable_malloc
!dec$ attributes alias:'_kmp_sharable_calloc_'::kmp_sharable_calloc
!dec$ attributes alias:'_kmp_sharable_realloc_'::kmp_sharable_realloc
!dec$ attributes alias:'_kmp_sharable_free_'::kmp_sharable_free
!dec$ attributes alias:'_kmp_lock_cond_wait_'::kmp_lock_cond_wait
!dec$ attributes alias:'_kmp_lock_cond_signal_'::kmp_lock_cond_signal
!dec$ attributes alias:'_kmp_lock_cond_broadcast_'::kmp_lock_cond_broadcast
!dec$ attributes alias:'_kmp_nest_lock_cond_wait_'::kmp_nest_lock_cond_wait
!dec$ attributes alias:'_kmp_nest_lock_cond_signal_'::kmp_nest_lock_cond_signal
!dec$ attributes alias:'_kmp_nest_lock_cond_broadcast_'::kmp_nest_lock_cond_broadcast
!dec$ attributes alias:'_kmp_get_num_processes_'::kmp_get_num_processes
!dec$ attributes alias:'_kmp_get_process_num_'::kmp_get_process_num
!dec$ attributes alias:'_kmp_get_process_thread_num_'::kmp_get_process_thread_num
!dec$ attributes alias:'_kmp_set_warnings_on_'::kmp_set_warnings_on
!dec$ attributes alias:'_kmp_set_warnings_off_'::kmp_set_warnings_off
!dec$ attributes alias:'_kmp_is_sharable_'::kmp_is_sharable
!dec$ attributes alias:'_kmp_deferred_atomic_add_i4_'::kmp_deferred_atomic_add_i4
!dec$ attributes alias:'_kmp_deferred_atomic_add_i8_'::kmp_deferred_atomic_add_i8
!dec$ attributes alias:'_kmp_deferred_atomic_add_r4_'::kmp_deferred_atomic_add_r4
!dec$ attributes alias:'_kmp_deferred_atomic_add_r8_'::kmp_deferred_atomic_add_r8

!dec$ endif

        end module omp_lib

