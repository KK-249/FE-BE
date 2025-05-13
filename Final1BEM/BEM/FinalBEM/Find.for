      subroutine Find_Element(NODE,m,n, NOD, J, K, found)
      integer NODE(m, n)
      integer NOD, J,K
      logical found
      integer ii, jj

      found = .false.
      do ii = 1, m
         do jj = 1, n
            if (NODE(ii,jj) .eq. NOD) then
               J = ii
               K = jj
               found = .true.
               return
            end if
         end do
      end do

      return
      end