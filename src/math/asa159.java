/*//
 * algorithm 159
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package math;

/**
 *
 * @author bukowski
 */
public class asa159 {
    
    public int seed;
    public int ierror;
    private int nrow;
    private int ncol;
    private int [] nrowt;
    private int [] ncolt;
    private int ntotal;
    private double [] fact;
    private int [] jwork;
            
    
    public asa159(int nrowIn, int ncolIn, int [] nrowtIn, int [] ncoltIn, int seedIn)
    {
        seed = seedIn;
        nrow = nrowIn;
        ncol = ncolIn;
        nrowt = nrowtIn;
        ncolt = ncoltIn;
        
        if ( nrow <= 1 )
        {
            //System.out.print ( "\n" );
            //System.out.print ( "RCONT - Fatal error!\n" );
            //System.out.print ( "  Input number of rows is less than 2.\n" );
            ierror = 1;
        }

    if ( ncol <= 1 )
    {
      //System.out.println ( "\n" );
      //System.out.print ( "RCONT - Fatal error!\n" );
      //System.out.print ( "  The number of columns is less than 2.\n" );
      ierror = 2;
    }

    for ( int i = 0; i < nrow; i++ )
    {
      if ( nrowt[i] <= 0 )
      {
        //System.out.print ( "\n" );
        //System.out.print ( "RCONT - Fatal error!\n" );
        //System.out.print ( "  An entry in the row sum vector is not positive.\n" );
        ierror = 3;
      }
    }

    for ( int j = 0; j < ncol; j++ )
    {
      if ( ncolt[j] <= 0 )
      {
        //System.out.print ( "\n" );
        //System.out.print ( "RCONT - Fatal error!\n" );
        //System.out.print ( "  An entry in the column sum vector is not positive.\n" );
        ierror = 4;
      }
    }

    if ( i4vec_sum ( ncol, ncolt ) != i4vec_sum ( nrow, nrowt ) )
    {
      //System.out.print ( "\n" );
      //System.out.print ( "RCONT - Fatal error!\n" );
      //System.out.print ( "  The row and column sum vectors do not have the same sum.\n" );
      ierror = 6;
    }

    ntotal = i4vec_sum ( ncol, ncolt );
    
    // log factorial array will need to be solved differently...
    // maybe in constructor
    fact = new double [ntotal + 1];
/*
  Calculate log-factorials.
*/
    double x = 0.0;
    fact[0] = 0.0;
    for ( int i = 1; i <= ntotal; i++ )
    {
      x = x + Math.log ( ( double ) ( i ) );
      fact[i] = x;
    }
    
    jwork = new int[ncol];

    } // end of constructor
    
    /******************************************************************************/

public double [][] rcont2 ()

/******************************************************************************/
/*
  Purpose:

    RCONT2 constructs a random two-way contingency table with given sums.

  Discussion:

    It is possible to specify row and column sum vectors which
    correspond to no table at all.  As far as I can see, this routine does
    not detect such a case.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 November 2010

  Author:

    Original FORTRAN77 version by WM Patefield.
    C version by John Burkardt.

  Reference:

    WM Patefield,
    Algorithm AS 159:
    An Efficient Method of Generating RXC Tables with
    Given Row and Column Totals,
    Applied Statistics,
    Volume 30, Number 1, 1981, pages 91-97.

  Parameters:

    Input, int NROW, NCOL, the number of rows and columns
    in the table.  NROW and NCOL must each be at least 2.

    Input, int NROWT[NROW], NCOLT[NCOL], the row and column
    sums.  Each entry must be positive.

    Input/output, int *KEY, a flag that indicates whether data has
    been initialized for this problem.  Set KEY = .FALSE. before the first
    call.

    Input/output, int *SEED, a seed for the random number generator.

    Output, int MATRIX[NROW*NCOL], the matrix.

    Output, int *IERROR, an error flag, which is returned
    as 0 if no error occurred.
*/
{
  int done1;
  int done2=0;
  int i;
  int ia;
  int iap;
  int ib=0;
  int ic;
  int id;
  int idp;
  int ie;
  int igp;
  int ihp;
  int ii;
  int iip;
  int j;
  int jc;
  int l;
  int lsm;
  int lsp;
  int m;
  int nll;
  int nlm;
  int nlmp;
  int nrowtl;
  double r;
  double sumprb;
  double x;
  double y;
  double [][] matrix = new double[nrow][ncol];
  boolean debug = false;

  ierror = 0;

/*
  Construct a random matrix.
*/
  for ( i = 0; i < ncol - 1; i++ )
  {
    jwork[i] = ncolt[i];
  }

  jc = ntotal;

  for ( l = 0; l < nrow - 1; l++ )
  {
    nrowtl = nrowt[l];
    ia = nrowtl;
    ic = jc;
    jc = jc - nrowtl;

    for ( m = 0; m < ncol - 1; m++ )
    {
        debug = false;
        if(m==73) { debug = true;}
        
      id = jwork[m];
      ie = ic;
      ic = ic - id;
      ib = ie - ia;
      ii = ib - id;
/*
  Test for zero entries in matrix.
*/
      if ( ie == 0 )
      {
        ia = 0;
        for ( j = m; j < ncol; j++ )
        {
          matrix[l][j] = 0;
        }
        break;
      }
/*
  Generate a pseudo-random number.
*/
      r = r8_uniform_01 ();
/*
  Compute the conditional expected value of MATRIX(L,M).
*/
      done1 = 0;

      for ( ; ; )
      {
        // sometimes ia*id goes over the limit of int32, resulting in negative numbers when cast onto double...
        //nlm = ( int ) ( ( double ) ( ia * id ) / ( double ) ( ie ) + 0.5 );
        nlm = ( int ) ( id * ( double ) ia / ( double ) ( ie ) + 0.5 );
        iap = ia + 1;
        idp = id + 1;
        igp = idp - nlm;
        ihp = iap - nlm;
        nlmp = nlm + 1;
        iip = ii + nlmp;
        x = Math.exp ( fact[iap-1] + fact[ib] + fact[ic] + fact[idp-1] -
          fact[ie] - fact[nlmp-1] - fact[igp-1] - fact[ihp-1] - fact[iip-1] );
        if ( r <= x )
        {
          break;
        }

        sumprb = x;
        y = x;
        nll = nlm;
        lsp = 0;
        lsm = 0;
/*
  Increment entry in row L, column M.
*/
        while ( lsp==0 )
        {
          j = ( id - nlm ) * ( ia - nlm );

          if ( j == 0 )
          {
            lsp = 1;
          }
          else
          {
            nlm = nlm + 1;
            x = x * ( double ) ( j ) / ( double ) ( nlm * ( ii + nlm ) );
            sumprb = sumprb + x;

            if ( r <= sumprb )
            {
              done1 = 1;
              break;
            }
          }

          done2 = 0;

          while ( lsm==0 )
          {
/*
  Decrement the entry in row L, column M.
*/
            j = nll * ( ii + nll );

            if ( j == 0 )
            {
              lsm = 1;
              break;
            }

            nll = nll - 1;
            y = y * ( double ) ( j ) / ( double ) ( ( id - nll ) * ( ia - nll ) );
            sumprb = sumprb + y;

            if ( r <= sumprb )
            {
              nlm = nll;
              done2 = 1;
              break;
            }

            if ( lsp==0 )
            {
              break;
            }

          }

          if ( done2==1 )
          {
            break;
          }

        }

        if ( done1==1 )
        {
          break;
        }

        if ( done2==1 )
        {
          break;
        }

        r = r8_uniform_01 ();
        r = sumprb * r;

      }

      matrix[l][m] = nlm;
      ia = ia - nlm;
      jwork[m] = jwork[m] - nlm;

    }
    matrix[l][(ncol-1)] = ia;
  }
/*
  Compute the last row.
*/
  for ( j = 0; j < ncol - 1; j++ )
  {
    matrix[nrow-1][j] = jwork[j];
  }
  matrix[nrow-1][(ncol-1)] = ib - matrix[nrow-1][(ncol-2)];

  return matrix;
}

private static int i4vec_sum(Integer n, int [] vec)
{
    int sum = 0;
    for(int i=0;i<n;i++)
    {
        sum += vec[i];
    }
    return sum;
}

private double r8_uniform_01()
{
    /******************************************************************************/
/*
  Purpose:

    R8_UNIFORM_01 returns a unit pseudorandom R8.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      r8_uniform_01 = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

    If the initial seed is 12345, then the first three computations are

      Input     Output      R8_UNIFORM_01
      SEED      SEED

         12345   207482415  0.096616
     207482415  1790989824  0.833995
    1790989824  2035175616  0.947702

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 August 2004

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation
    edited by Jerry Banks,
    Wiley Interscience, page 95, 1998.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.

    P A Lewis, A S Goodman, J M Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.

  Parameters:

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has been updated.

    Output, double R8_UNIFORM_01, a new pseudorandom variate, strictly between
    0 and 1.
*/

  int k;
  double r;

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + 2147483647;
  }
/*
  Although SEED can be represented exactly as a 32 bit integer,
  it generally cannot be represented exactly as a 32 bit real number!
*/
  r = ( ( double ) ( seed ) ) * 4.656612875E-10;

  return r;
}
    
}
