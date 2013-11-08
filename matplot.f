C
C
        SUBROUTINE RSPLOT(X,Y,N,NBOD,IW)
        IMPLICIT REAL *8 (A-H,O-Z)
        DIMENSION X(1),Y(1),V(4)
C
C        THIS SUBROUTINE PLOTS THE CURVE SPECIFIED BY ITS NODES
C        X,Y (CONCEPTUALLY, THE CURVE CONNECTS THE POINT (X(1),Y(1))
C        WITH THE POINT (X(N),Y(N)). 
C
C  INPUT PARAMETERS:
C   X,Y - THE COORDINATES OF THE CURVE TO PLOT
C   N - THE NUMBER OF ELEMENTS IN ARAYS X,Y.
C   IW - THE FORTRAN UNIT NUMBER ON WHICH THE OUTPUT DATA SET IS WRITTEN
C
C  OUTPUT PARAMETERS: NONE
C
C
      WRITE(IW,*) 'x = ['
      WRITE(IW,1400) (X(I),Y(I),I=1,N)
ccc      write (iw,1400) x(1),y(1)
      WRITE(IW,*) '];'
1300  FORMAT(D15.6)
1400  FORMAT(2D15.6)
      write(iw,*) ' xx = x(:,1);'
      write(iw,*) ' yy = x(:,2);'
      write(iw, *) ' plot(xx,yy,''k'')'
ccc      write(iw, *) ' plot(xx,yy," o \')'
      write(iw,*) 'axis equal'
ccc      write(iw, *) ' plot(x)'
      IF (NBOD.eq.1) write(iw,*) ' hold on'
      RETURN
C
      entry rsplt1(x,n,iw)
      WRITE(IW,*) 'x = ['
      WRITE(IW,1300) (X(I),I=1,N)
      WRITE(IW,*) '];'
      RETURN
C
C
      ENTRY RSPINI(IW,XMIN,XMAX,YMIN,YMAX)
      II = 0
      V(1) = XMIN
      V(2) = XMAX
      V(3) = YMIN
      V(4) = YMAX
      write(IW,*) 'v = ['
      write(iw,1300) (v(i),i=1,4)
      write(iw,*) '];'
      write(iw,*) ' axis(v)'
ccc      write(iw,*) ' axis(\'square\')'
          RETURN
          END
c
        SUBROUTINE RS_3D_PLOT(X,Y,Z,N,NBOD,IW)
        IMPLICIT REAL *8 (A-H,O-Z)
        DIMENSION X(1),Y(1),z(1), V(4)
C
C        THIS SUBROUTINE PLOTS THE CURVE SPECIFIED BY ITS NODES
C        X,Y (CONCEPTUALLY, THE CURVE CONNECTS THE POINT (X(1),Y(1))
C        WITH THE POINT (X(N),Y(N)). 
C
C  INPUT PARAMETERS:
C   X,Y - THE COORDINATES OF THE CURVE TO PLOT
C   N - THE NUMBER OF ELEMENTS IN ARAYS X,Y.
C   IW - THE FORTRAN UNIT NUMBER ON WHICH THE OUTPUT DATA SET IS WRITTEN
C
C  OUTPUT PARAMETERS: NONE
C
C
      WRITE(IW,*) 'x = ['
      WRITE(IW,1400) (X(I),Y(I),Z(I),I=1,N)
      write (iw,1400) x(1),y(1),z(1)
      WRITE(IW,*) '];'
1300  FORMAT(D15.6)
1400  FORMAT(3D15.6)
      write(iw,*) ' xx = x(:,1);'
      write(iw,*) ' yy = x(:,2);'
      write(iw,*) ' zz = x(:,3);'
ccc      write(iw, *) ' plot(xx,yy)'
      write(iw, *) ' plot3(xx,yy,zz,''k'',''LineWidth'',2)'
ccc      write(iw,*) 'axis equal'
ccc      write(iw, *) ' plot(x)'
      IF (NBOD.eq.1) write(iw,*) ' hold on'
      RETURN
      end
C
C
        SUBROUTINE RSLOGPLOT(X,Y,N,NBOD,IW)
        IMPLICIT REAL *8 (A-H,O-Z)
        DIMENSION X(1),Y(1),V(4)
C
C        THIS SUBROUTINE PLOTS THE CURVE SPECIFIED BY ITS NODES
C        X,Y (CONCEPTUALLY, THE CURVE CONNECTS THE POINT (X(1),Y(1))
C        WITH THE POINT (X(N),Y(N)). 
C
C  INPUT PARAMETERS:
C   X,Y - THE COORDINATES OF THE CURVE TO PLOT
C   N - THE NUMBER OF ELEMENTS IN ARAYS X,Y.
C   IW - THE FORTRAN UNIT NUMBER ON WHICH THE OUTPUT DATA SET IS WRITTEN
C
C  OUTPUT PARAMETERS: NONE
C
C
      WRITE(IW,*) 'x = ['
      WRITE(IW,1400) (X(I),Y(I),I=1,N)
ccc      write (iw,1400) x(1),y(1)
      WRITE(IW,*) '];'
1300  FORMAT(D15.6)
1400  FORMAT(2D15.6)
      write(iw,*) ' xx = x(:,1);'
      write(iw,*) ' yy = x(:,2);'
      write(iw, *) ' semilogy(xx,yy)'
ccc      write(iw, *) ' plot(xx,yy,\'o\')'
ccc      write(iw,*) 'axis equal'
ccc      write(iw, *) ' plot(x)'
      IF (NBOD.eq.1) write(iw,*) ' hold on'
      RETURN
      end
C
C
        SUBROUTINE RSCPLOT(Z,N,NBOD,IW)
        IMPLICIT REAL *8 (A-H,O-Z)
        DIMENSION V(4)
        complex*16 z(n)
C
C        THIS SUBROUTINE PLOTS THE CURVE SPECIFIED BY ITS NODES
C        z = x + Iy (CONCEPTUALLY, THE CURVE CONNECTS THE POINT (X(1),Y(1))
C        WITH THE POINT (X(N),Y(N)). 
C
C  INPUT PARAMETERS:
C   X,Y - THE COORDINATES OF THE CURVE TO PLOT
C   N - THE NUMBER OF ELEMENTS IN ARAYS X,Y.
C   IW - THE FORTRAN UNIT NUMBER ON WHICH THE OUTPUT DATA SET IS WRITTEN
C
C  OUTPUT PARAMETERS: NONE
C
C
      WRITE(IW,*) 'x = ['
      WRITE(IW,1400) (z(I),I=1,N)
ccc      write (iw,1400) x(1),y(1)
      WRITE(IW,*) '];'
1300  FORMAT(D15.6)
1400  FORMAT(2D15.6)
      write(iw,*) ' xx = x(:,1);'
      write(iw,*) ' yy = x(:,2);'
      if (nbod.eq.1) then
      write(iw, *) ' plot(xx,yy,''k'',''LineWidth'',2)'
       else
      write(iw, *) ' fill(xx,yy,''k'')'       
      end if
ccc      write(iw, *) ' plot(xx,yy,'' . '' )'
      write(iw,*) 'axis equal'
ccc      write(iw, *) ' plot(x)'
      IF (NBOD.eq.1) write(iw,*) ' hold on'
      RETURN
      end
c
        SUBROUTINE RSCPLOT_DOT(Z,N,NBOD,IW)
        IMPLICIT REAL *8 (A-H,O-Z)
        DIMENSION V(4)
        complex*16 z(n)
C
C        THIS SUBROUTINE PLOTS THE CURVE SPECIFIED BY ITS NODES
C        z = x + Iy (CONCEPTUALLY, THE CURVE CONNECTS THE POINT (X(1),Y(1))
C        WITH THE POINT (X(N),Y(N)). 
C
C  INPUT PARAMETERS:
C   X,Y - THE COORDINATES OF THE CURVE TO PLOT
C   N - THE NUMBER OF ELEMENTS IN ARAYS X,Y.
C   IW - THE FORTRAN UNIT NUMBER ON WHICH THE OUTPUT DATA SET IS WRITTEN
C
C  OUTPUT PARAMETERS: NONE
C
C
      WRITE(IW,*) 'x = ['
      WRITE(IW,1400) (z(I),I=1,N)
ccc      write (iw,1400) x(1),y(1)
      WRITE(IW,*) '];'
1300  FORMAT(D15.6)
1400  FORMAT(2D15.6)
      write(iw,*) ' xx = x(:,1);'
      write(iw,*) ' yy = x(:,2);'
      write(iw, *) ' plot(xx,yy,''r.'',''MarkerSize'',0.5)'
ccc      write(iw, *) ' plot(xx,yy,'' . '' )'
      write(iw,*) 'axis equal'
ccc      write(iw, *) ' plot(x)'
      IF (NBOD.eq.1) write(iw,*) ' hold on'
      RETURN
      end
c
        SUBROUTINE RSC_STAR_PLOT(Z,N,IW)
        IMPLICIT REAL *8 (A-H,O-Z)
        DIMENSION V(4)
        complex*16 z(n)
C
C        THIS SUBROUTINE PLOTS THE CURVE SPECIFIED BY ITS NODES
C        z = x + Iy (CONCEPTUALLY, THE CURVE CONNECTS THE POINT (X(1),Y(1))
C        WITH THE POINT (X(N),Y(N)). 
C
C  INPUT PARAMETERS:
C   X,Y - THE COORDINATES OF THE CURVE TO PLOT
C   N - THE NUMBER OF ELEMENTS IN ARAYS X,Y.
C   IW - THE FORTRAN UNIT NUMBER ON WHICH THE OUTPUT DATA SET IS WRITTEN
C
C  OUTPUT PARAMETERS: NONE
C
C
      WRITE(IW,*) 'x = ['
      WRITE(IW,1400) (z(I),I=1,N)
ccc      write (iw,1400) x(1),y(1)
      WRITE(IW,*) '];'
1300  FORMAT(D15.6)
1400  FORMAT(2D15.6)
      write(iw,*) ' xx = x(:,1);'
      write(iw,*) ' yy = x(:,2);'
      write(iw, *) ' plot(xx,yy,''r*'')'
ccc      write(iw, *) ' plot(xx,yy,'' . '' )'
      write(iw,*) 'axis equal'
ccc      write(iw, *) ' plot(x)'
      IF (NBOD.eq.1) write(iw,*) ' hold on'
      RETURN
      end
      
