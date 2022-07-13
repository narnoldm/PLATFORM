#include "extern_func.hpp"

using namespace :: std;

// need these because C++ treats double quotes as string literals
char ROW[1] = {'R'};
char COLUMN[1] = {'C'};
char UNPACKING[1] = {'U'};
char PACKING[1] = {'P'};
char TRAN[1] = {'T'};
char NOTRAN[1] = {'N'};
char TOP_GET[1] = {'!'};
char BCAST[1] = {'B'};

// multiply A and B, scaled by ALPHA, store result in B
void dmmmult_(int M, int N, char* ALPHA, char* A, int LDA, char* B, int LDB)
{

    double one, zero;
    one = 1.0e+0;
    zero = 0.0e+0;

    double alpha = *((double*)ALPHA);
    double* a = (double*)A;
    double* b = (double*)B;

    int idxa, idxb;

    if (alpha == one)
    {
        for (int j = 0; j < N; ++j)
        {
            for (int i = 0; i < M; ++i)
            {
                idxa = j * LDA + i;
                idxb = j * LDB + i;
                b[idxb] = a[idxa] * b[idxb];
            }
        }
    }
    else
    {
        for (int j = 0; j < N; ++j)
        {
            for (int i = 0; i < M; ++i)
            {
                idxa = j * LDA + i;
                idxb = j * LDB + i;
                b[idxb] = alpha * a[idxa] * b[idxb];
            }
        }
    }
}

// multiply A and B, scaled by ALPHA, store result in A
void dmmtlum_(int M, int N, char* ALPHA, char* A, int LDA, char* B, int LDB)
{

    double one, zero;
    one = 1.0e+0;
    zero = 0.0e+0;

    double alpha = *((double*)ALPHA);
    double* a = (double*)A;
    double* b = (double*)B;

    int idxa, idxb;

    if (alpha == one)
    {
        for (int j = 0; j < N; ++j)
        {
            for (int i = 0; i < M; ++i)
            {
                idxa = j * LDA + i;
                idxb = j * LDB + i;
                a[idxb] = a[idxa] * b[idxb];
            }
        }
    }
    else
    {
        for (int j = 0; j < N; ++j)
        {
            for (int i = 0; i < M; ++i)
            {
                idxa = j * LDA + i;
                idxb = j * LDB + i;
                a[idxb] = alpha * a[idxa] * b[idxb];
            }
        }
    }
}


void pdgemul(char trans, int m, int n,
            double alpha, double* AMat, int ai, int aj, int* adesc,
            double* CMat, int ci, int cj, int* cdesc)
{

    char        ACroc, * one, * talpha, * zero;
    int         ACmyprocD, ACmyprocR, ACnD, ACnR, ACnprocsD, ACnprocsR,
                Abufld, AcurrocR, Afr, Afwd, AiD, AiR, AiiD, AiiR, AinbD,
                AinbR, Ainb1D, Ainb1R, AisR, Akk, Ald, AnbD, AnbR, AnpD,
                AnpR, Aoff, ArocD, ArocR, AsrcR, Cbufld, CcurrocR, Cfr,
                Cfwd, CiD, CiR, CiiD, CiiR, CinbD, CinbR, Cinb1D, Cinb1R,
                CisR, Ckk, Cld, CnbD, CnbR, CnpD, CnpR, Coff, CrocD, CrocR,
                CsrcR, ctxt, k, kb, kbb, lcmb, maxp, maxpm1, maxpq, maxq,
                mycol, myrow, npcol, npq, nprow, ncpq, nrpq, p=0, q=0,
                row2row, size, tmp;
    PB_VM_T     VM;
    PBTYP_T* TYPE = PB_Cdtypeset();

    /*
    *  .. Local Arrays ..
    */
    int DBUFA[DLEN_], DBUFC[DLEN_];
    char * Abuf = NULL, * Cbuf = NULL;
    char ALPHA[sizeof(alpha)], BETA[sizeof(double)];
    memcpy(ALPHA, &alpha, sizeof(alpha));
    double beta = 0.0;
    memcpy(BETA, &beta, sizeof(double));

    double *A_arr, *C_arr;
    char* A = (char *) AMat;
    char* C = (char *) CMat;

    /*
    *  Retrieve process grid information
    */
    int ia, ja, ic, jc;
    int desca[DLEN_], descc[DLEN_];
    PB_CargFtoC( ai, aj, adesc, &ia, &ja, desca );
    PB_CargFtoC( ci, cj, cdesc, &ic, &jc, descc );
    Cblacs_gridinfo( ( ctxt = descc[CTXT_] ), &nprow, &npcol, &myrow, &mycol );

    /*
    *  Loop over the rows of sub( C ) when m <= n, and the columns of sub( C )
    *  otherwise.
    */
    row2row = ( ( m <= n ) || ( npcol == 1 ) || ( desca[CSRC_] == -1 ) );
    if( row2row )
    {
        AinbR = desca[IMB_]; AnbR = desca[MB_]; AsrcR = desca[RSRC_];
        CinbR = descc[IMB_]; CnbR = descc[MB_]; CsrcR = descc[RSRC_];

        /*
        *  If sub( A ) and sub( C ) span only one process row, then there is no need
        *  to pack the data.
        */
        if( !( PB_Cspan( m, ia, AinbR, AnbR, AsrcR, nprow ) ) &&
            !( PB_Cspan( m, ic, CinbR, CnbR, CsrcR, nprow ) ) )
        {
            daxy( TYPE, &trans, m, n, ALPHA, A, ia, ja, desca, ROW,
                        C, ic, jc, descc, ROW);
            return;
        }

        /*
        *  Compute local information for sub( A ) and sub( C )
        */
        ACnR      = m;           ACnD      = n;
        ACmyprocR = myrow;       ACnprocsR = nprow;
        ACmyprocD = mycol;       ACnprocsD = npcol;      ACroc = 'R';
        AiR       = ia;          AiD       = ja;
        AinbD     = desca[INB_]; AnbD      = desca[NB_]; Ald   = desca[LLD_];
        PB_Cinfog2l( ia, ja, desca, ACnprocsR, ACnprocsD, ACmyprocR, ACmyprocD,
                        &AiiR, &AiiD, &ArocR, &ArocD );
        CiR       = ic;          CiD       = jc;
        CinbD     = descc[INB_]; CnbD      = descc[NB_]; Cld = descc[LLD_];
        PB_Cinfog2l( ic, jc, descc, ACnprocsR, ACnprocsD, ACmyprocR, ACmyprocD,
                    &CiiR, &CiiD, &CrocR, &CrocD );
    }
    else
    {
        AinbR = desca[INB_]; AnbR = desca[NB_]; AsrcR = desca[CSRC_];
        CinbR = descc[INB_]; CnbR = descc[NB_]; CsrcR = descc[CSRC_];

        /*
        *  If sub( A ) and sub( C ) span only one process column, then there is no need
        *  to pack the data.
        */
        if( !( PB_Cspan( n, ja, AinbR, AnbR, AsrcR, npcol ) ) &&
            !( PB_Cspan( n, jc, CinbR, CnbR, CsrcR, npcol ) ) )
        {
            daxy( TYPE, &trans, m, n, ALPHA, A, ia, ja, desca, COLUMN,
                        C, ic, jc, descc, COLUMN );
            return;
        }

        /*
        *  Compute local information for sub( A ) and sub( C )
        */
        ACnR      = n;           ACnD      = m;
        ACmyprocR = mycol;       ACnprocsR = npcol;
        ACmyprocD = myrow;       ACnprocsD = nprow;      ACroc = 'C';
        AiR       = ja;          AiD       = ia;
        AinbD     = desca[IMB_]; AnbD      = desca[MB_]; Ald   = desca[LLD_];
        PB_Cinfog2l( ia, ja, desca, ACnprocsD, ACnprocsR, ACmyprocD, ACmyprocR,
                        &AiiD, &AiiR, &ArocD, &ArocR );
        CiR       = jc;          CiD       = ic;
        CinbD     = descc[IMB_]; CnbD      = descc[MB_]; Cld = descc[LLD_];
        PB_Cinfog2l( ic, jc, descc, ACnprocsD, ACnprocsR, ACmyprocD, ACmyprocR,
                    &CiiD, &CiiR, &CrocD, &CrocR );
    }

    size   = TYPE->size; one = TYPE->one; zero = TYPE->zero;
    kb     = pilaenv_( &ctxt, C2F_CHAR( &TYPE->type ) );

    Ainb1D = PB_Cfirstnb( ACnD, AiD, AinbD, AnbD );
    AnpD   = PB_Cnumroc( ACnD, 0, Ainb1D, AnbD, ACmyprocD, ArocD, ACnprocsD );
    Ainb1R = PB_Cfirstnb( ACnR, AiR, AinbR, AnbR );
    AisR   = ( ( AsrcR < 0 ) || ( ACnprocsR == 1 ) );

    Cinb1D = PB_Cfirstnb( ACnD, CiD, CinbD, CnbD );
    CnpD   = PB_Cnumroc( ACnD, 0, Cinb1D, CnbD, ACmyprocD, CrocD, ACnprocsD );
    Cinb1R = PB_Cfirstnb( ACnR, CiR, CinbR, CnbR );
    CisR   = ( ( CsrcR < 0 ) || ( ACnprocsR == 1 ) );

    lcmb   = PB_Clcm( ( maxp = ( CisR ? 1 : ACnprocsR ) ) * CnbR,
                      ( maxq = ( AisR ? 1 : ACnprocsR ) ) * AnbR );

    Afwd   = 1;
    Cfwd   = 1;

    /*
    *  When sub( A ) is not replicated and backward pass on sub( A ), find the
    *  virtual process q owning the last row or column of sub( A ).
    */
    if( !( AisR ) && !( Afwd ) )
    {
        tmp = PB_Cindxg2p( ACnR-1, Ainb1R, AnbR, ArocR, ArocR, ACnprocsR );
        q   = MModSub( tmp, ArocR, ACnprocsR );
    }

    /*
    *  When sub( C ) is not replicated and backward pass on sub( C ), find the
    *  virtual process p owning the last row or column of sub( C ).
    */
    if( !( CisR ) && !( Cfwd ) )
    {
        tmp = PB_Cindxg2p( ACnR-1, Cinb1R, CnbR, CrocR, CrocR, ACnprocsR );
        p   = MModSub( tmp, CrocR, ACnprocsR );
    }

    /*
    *  Loop over the processes of the virtual grid
    */
    maxpm1 = maxp - 1; maxpq = maxp * maxq;
    for( k = 0; k < maxpq; k++ )
    {
        AcurrocR = ( AisR ? -1 : MModAdd( ArocR, q, ACnprocsR ) );
        CcurrocR = ( CisR ? -1 : MModAdd( CrocR, p, ACnprocsR ) );

        if( ( AisR || ( ACmyprocR == AcurrocR ) ) ||
            ( CisR || ( ACmyprocR == CcurrocR ) ) )
        {
            Ckk = CiiR; Akk = AiiR;

            /*
            *  Initialize local virtual matrix in process (p,q)
            */
            AnpR = PB_Cnumroc( ACnR, 0, Ainb1R, AnbR, AcurrocR, ArocR, ACnprocsR );
            CnpR = PB_Cnumroc( ACnR, 0, Cinb1R, CnbR, CcurrocR, CrocR, ACnprocsR );
            PB_CVMinit( &VM, 0, CnpR, AnpR, Cinb1R, Ainb1R, CnbR, AnbR, p, q,
                        maxp, maxq, lcmb );

            /*
            *  Figure out how many diagonal entries in this new virtual process (npq).
            */
            npq = PB_CVMnpq( &VM );

            /*
            *  Re-adjust the number of rows or columns to be (un)packed, in order to average
            *  the message sizes.
            */
            if( npq ) kbb = npq / ( ( npq - 1 ) / kb + 1 );

            if( row2row )
            {
                while( npq )
                {
                    kbb = min( kbb, npq );

                    /*
                    *  Find out how many rows of sub( A ) and sub( C ) are contiguous
                    */
                    PB_CVMcontig( &VM, &nrpq, &ncpq, &Coff, &Aoff );

                    /*
                    *  Compute the descriptor DBUFA for the buffer that will contained the packed
                    *  rows of sub( A ).
                    */
                    if( ( Afr = ( ncpq < kbb ) ) != 0 )
                    {
                        /*
                        *  If rows of sub( A ) are not contiguous, then allocate the buffer and pack
                        *  the kbb rows of sub( A ).
                        */
                        Abufld = kbb;
                        if( AisR || ( ACmyprocR == AcurrocR ) )
                        {
                            Abuf = PB_Cmalloc( AnpD * kbb * size );
                            PB_CVMpack( TYPE, &VM, COLUMN, &ACroc, PACKING, NOTRAN,
                                        kbb, AnpD, one, Mptr( A, Akk, AiiD, Ald,
                                        size ), Ald, zero,  Abuf, Abufld );
                        }
                    }
                    else
                    {
                        /*
                        *  Otherwise, re-use sub( A ) directly.
                        */
                        Abufld = Ald;
                        if( AisR || ( ACmyprocR == AcurrocR ) )
                            Abuf = Mptr( A, Akk+Aoff, AiiD, Ald, size );
                    }
                    PB_Cdescset( DBUFA, kbb, ACnD, kbb, Ainb1D, kbb, AnbD, AcurrocR,
                                ArocD, ctxt, Abufld );

                    /*
                    *  Compute the descriptor DBUFC for the buffer that will contained the packed
                    *  rows of sub( C ). Allocate it.
                    */
                    if( ( Cfr = ( nrpq < kbb ) ) != 0 )
                    {
                        /*
                        *  If rows of sub( C ) are not contiguous, then allocate receiving buffer.
                        */
                        Cbufld = kbb; talpha = one;
                        if( CisR || ( ACmyprocR == CcurrocR ) )
                            Cbuf = PB_Cmalloc( CnpD * kbb * size );
                    }
                    else
                    {
                        /*
                        *  Otherwise, re-use sub( C ) directly.
                        */
                        Cbufld = Cld; talpha = ALPHA;
                        if( CisR || ( ACmyprocR == CcurrocR ) )
                            Cbuf = Mptr( C, Ckk+Coff, CiiD, Cld, size );
                    }
                    PB_Cdescset( DBUFC, kbb, ACnD, kbb, Cinb1D, kbb, CnbD, CcurrocR,
                                CrocD, ctxt, Cbufld );

                    /*
                    *  Add the one-dimensional buffer Abuf into Cbuf.
                    */
                    daxy( TYPE, &trans, kbb, ACnD, talpha, Abuf, 0, 0, DBUFA,
                                &ACroc, Cbuf, 0, 0, DBUFC, &ACroc );

                    /*
                    *  Release the buffer containing the packed rows of sub( A )
                    */
                    if( Afr && ( AisR || ( ACmyprocR == AcurrocR ) ) )
                        if( Abuf ) free( Abuf );

                    /*
                    *  Unpack the kbb rows of sub( C ) and release the buffer containing them.
                    */
                    if( Cfr && ( CisR || ( ACmyprocR == CcurrocR ) ) )
                    {
                        // TODO: does this need to be converted?
                        cout << "ERROR: Chris didn't implement this" << endl;
                        MPI_Abort(MPI_COMM_WORLD, -1);
                        PB_CVMpack( TYPE, &VM, ROW, &ACroc, UNPACKING, NOTRAN, kbb,
                                    CnpD, BETA, Mptr( C, Ckk, CiiD, Cld, size ), Cld,
                                    ALPHA, Cbuf, Cbufld );
                        if( Cbuf ) free( Cbuf );
                    }

                    /*
                    *  Update the local row indexes of sub( A ) and sub( C )
                    */
                    PB_CVMupdate( &VM, kbb, &Ckk, &Akk );
                    npq -= kbb;
                }
            }
            else
            {
                while( npq )
                {
                    kbb = min( kbb, npq );

                    /*
                    *  Find out how many columns of sub( A ) and sub( C ) are contiguous
                    */
                    PB_CVMcontig( &VM, &nrpq, &ncpq, &Coff, &Aoff );

                    /*
                    *  Compute the descriptor DBUFA for the buffer that will contained the packed
                    *  columns of sub( A ).
                    */
                    if( ( Afr = ( ncpq < kbb ) ) != 0 )
                    {
                        /*
                        *  If columns of sub( A ) are not contiguous, then allocate the buffer and
                        *  pack the kbb columns of sub( A ).
                        */
                        Abufld = max( 1, AnpD );
                        if( AisR || ( ACmyprocR == AcurrocR ) )
                        {
                            Abuf = PB_Cmalloc( AnpD * kbb * size );
                            PB_CVMpack( TYPE, &VM, COLUMN, &ACroc, PACKING, NOTRAN,
                                        kbb, AnpD, one, Mptr( A, AiiD, Akk, Ald,
                                        size ), Ald, zero,  Abuf, Abufld );
                        }
                    }
                    else
                    {
                        /*
                        *  Otherwise, re-use sub( A ) directly.
                        */
                        Abufld = Ald;
                        if( AisR || ( ACmyprocR == AcurrocR ) )
                            Abuf = Mptr( A, AiiD, Akk+Aoff, Ald, size );
                    }
                    PB_Cdescset( DBUFA, ACnD, kbb, Ainb1D, kbb, AnbD, kbb, ArocD,
                                AcurrocR, ctxt, Abufld );

                    /*
                    *  Compute the descriptor DBUFC for the buffer that will contained the packed
                    *  columns of sub( C ). Allocate it.
                    */
                    if( ( Cfr = ( nrpq < kbb ) ) != 0 )
                    {
                        /*
                        *  If columns of sub( C ) are not contiguous, then allocate receiving buffer.
                        */
                        Cbufld = max( 1, CnpD ); talpha = one;
                        if( CisR || ( ACmyprocR == CcurrocR ) )
                            Cbuf = PB_Cmalloc( CnpD * kbb * size );
                    }
                    else
                    {
                        Cbufld = Cld; talpha = ALPHA;
                        if( CisR || ( ACmyprocR == CcurrocR ) )
                            Cbuf = Mptr( C, CiiD, Ckk+Coff, Cld, size );
                    }
                    PB_Cdescset( DBUFC, ACnD, kbb, Cinb1D, kbb, CnbD, kbb, CrocD,
                                CcurrocR, ctxt, Cbufld );

                    /*
                    *  Add the one-dimensional buffer Abuf into Cbuf.
                    */
                    daxy( TYPE, &trans, ACnD, kbb, talpha, Abuf, 0, 0, DBUFA,
                                &ACroc, Cbuf, 0, 0, DBUFC, &ACroc );

                    /*
                    *  Release the buffer containing the packed columns of sub( A )
                    */
                    if( Afr && ( AisR || ( ACmyprocR == AcurrocR ) ) )
                        if( Abuf ) free( Abuf );

                    /*
                    *  Unpack the kbb columns of sub( C ) and release the buffer containing them.
                    */
                    if( Cfr && ( CisR || ( ACmyprocR == CcurrocR ) ) )
                    {
                        // TODO: does this need to be converted?
                        cout << "ERROR: Chris didn't implement this" << endl;
                        MPI_Abort(MPI_COMM_WORLD, -1);
                        PB_CVMpack( TYPE, &VM, ROW, &ACroc, UNPACKING, NOTRAN, kbb,
                                    CnpD, BETA, Mptr( C, CiiD, Ckk, Cld, size ), Cld,
                                    ALPHA, Cbuf, Cbufld );
                        if( Cbuf ) free( Cbuf );
                    }

                    /*
                    *  Update the local row index of sub( A ) and the local column index of sub( C )
                    */
                    PB_CVMupdate( &VM, kbb, &Ckk, &Akk );
                    npq -= kbb;
                }
            }
        }

        /*
        *  Go to the next virtual process (p,q)
        */
        if( ( Cfwd && ( p == maxpm1 ) ) || ( !( Cfwd ) && ( p == 0 ) ) )
            q = ( Afwd ? MModAdd1( q, maxq ) : MModSub1( q, maxq ) );
        p = ( Cfwd ? MModAdd1( p, maxp ) : MModSub1( p, maxp ) );
    }

}

void daxy( PBTYP_T* TYPE, char* CONJUG, int M, int N, char* ALPHA,
    char* A, int IA, int JA, int* DESCA, char* AROC,
    char* B, int IB, int JB, int* DESCB, char* BROC)
{

    /*
    *  .. Local Scalars ..
    */
    char    ascope, bscope, * buf = NULL, * one, * top, tran, * zero;
    int     Acol, Aii, AinbD, Ainb1D, AisD, AisR, AisRow, AiD, Ajj, Ald,
            AmyprocD, AmyprocR, AnbD, AnD, AnR, AnpD, AnprocsD, AnprocsR,
            AprocD, AprocR, Aroc, Arow, Bcol, Bii, BinbD, Binb1D, BisD,
            BisR, BisRow, BiD, Bjj, Bld, BmyprocD, BmyprocR, BnbD, BnD,
            BnR, BnpD, BnprocsD, BnprocsR, BprocD, BprocR, Broc, Brow,
            BsrcD, OneBlock, OneDgrid, RRorCC, Square, cdst, csrc, ctxt,
            dst, gcdPQ, k, l, lcmPQ, lcmb, ma, mb, mycol, myrow, na, nb,
            npcol, npq, nprow, p, q, rdst, rsrc, size, src;
    PB_VM_T VM;
    MMMULT_T mult;

    // TODO: generalize
    MMMULT_T dmmmult = dmmmult_;
    MMMULT_T dmmtlum = dmmtlum_;
    double beta = 0.0;
    char BETA[sizeof(double)];
    memcpy(BETA, &beta, sizeof(double));

    /*
    *  Quick return if possible
    */
    if( ( M <= 0 ) || ( N <= 0 ) ) return;

    /*
    *  Retrieve process grid information
    */
    Cblacs_gridinfo( ( ctxt = DESCA[ CTXT_ ] ), &nprow, &npcol, &myrow, &mycol );

    /*
    *  Determine if sub( A ) is distributed or not
    */
    if( ( AisRow = ( Mupcase( AROC[0] ) == CROW ) ) != 0 )
        AisD = ( ( DESCA[CSRC_] >= 0 ) && ( ( AnprocsD = npcol ) > 1 ) );
    else
        AisD = ( ( DESCA[RSRC_] >= 0 ) && ( ( AnprocsD = nprow ) > 1 ) );

    /*
    *  Determine if sub( B ) is distributed or not
    */
    if( ( BisRow = ( Mupcase( BROC[0] ) == CROW ) ) != 0 )
       BisD = ( ( DESCB[CSRC_] >= 0 ) && ( ( BnprocsD = npcol ) > 1 ) );
    else
       BisD = ( ( DESCB[RSRC_] >= 0 ) && ( ( BnprocsD = nprow ) > 1 ) );

    /*
    *  AisD && BisD <=> both operands are indeed distributed
    */
    if( AisD && BisD )
    {
        /*
        *  Retrieve sub( A )'s local information: Aii, Ajj, Arow, Acol ...
        */
        PB_Cinfog2l( IA, JA, DESCA, nprow, npcol, myrow, mycol, &Aii, &Ajj, &Arow,
                        &Acol );
        if( AisRow )
        {
            AinbD  = DESCA[INB_]; AnbD = DESCA[NB_]; Ald = DESCA[LLD_];
            AiD    = JA;    AnD    =  N;      AnR    =  M;
            AprocD = Acol; AmyprocD = mycol;
            AprocR = Arow; AmyprocR = myrow; AnprocsR = nprow;
            AisR   = ( ( DESCA[ RSRC_ ] == -1 ) || ( AnprocsR == 1 ) );
        }
        else
        {
            AinbD  = DESCA[IMB_]; AnbD = DESCA[MB_]; Ald = DESCA[LLD_];
            AiD    = IA;   AnD    =  M;      AnR    =  N;
            AprocD = Arow; AmyprocD = myrow;
            AprocR = Acol; AmyprocR = mycol; AnprocsR = npcol;
            AisR   = ( ( DESCA[ CSRC_ ] == -1 ) || ( AnprocsR == 1 ) );
        }
        Ainb1D = PB_Cfirstnb( AnD, AiD, AinbD, AnbD );

        /*
        *  Retrieve sub( B )'s local information: Bii, Bjj, Brow, Bcol ...
        */
        PB_Cinfog2l( IB, JB, DESCB, nprow, npcol, myrow, mycol, &Bii, &Bjj, &Brow,
                        &Bcol );
        if( BisRow )
        {
            BinbD  = DESCB[ INB_  ]; BnbD   = DESCB[ NB_   ];
            BsrcD  = DESCB[ CSRC_ ]; Bld    = DESCB[ LLD_  ];
            BiD    = JB;
            if( AisRow ) { BnD = N; BnR = M; } else { BnD = M; BnR = N; }
            BprocD = Bcol; BmyprocD = mycol;
            BprocR = Brow; BmyprocR = myrow; BnprocsR = nprow;
            BisR   = ( ( DESCB[ RSRC_ ] == -1 ) || ( BnprocsR == 1 ) );
        }
        else
        {
            BinbD  = DESCB[ IMB_  ]; BnbD   = DESCB[ MB_   ];
            BsrcD  = DESCB[ RSRC_ ]; Bld    = DESCB[ LLD_  ];
            BiD    = IB;
            if( AisRow ) { BnD = N; BnR = M; } else { BnD = M; BnR = N; }
            BprocD = Brow; BmyprocD = myrow;
            BprocR = Bcol; BmyprocR = mycol; BnprocsR = npcol;
            BisR   = ( ( DESCB[ CSRC_ ] == -1 ) || ( BnprocsR == 1 ) );
        }
        Binb1D = PB_Cfirstnb( BnD, BiD, BinbD, BnbD );

        /*
        *  Are sub( A ) and sub( B ) both row or column vectors ?
        */
        RRorCC = ( ( AisRow && BisRow ) || ( !( AisRow ) && !( BisRow ) ) );

        /*
        *  Do sub( A ) and sub( B ) span more than one process ?
        */
        OneDgrid = ( ( AnprocsD ==   1 ) && ( BnprocsD ==   1 ) );
        OneBlock = ( ( Ainb1D   >= AnD ) && ( Binb1D   >= BnD ) );

        /*
        *  Are sub( A ) and sub( B ) distributed in the same manner ?
        */
        Square   = ( ( Ainb1D   ==   Binb1D ) && ( AnbD == BnbD ) &&
                     ( AnprocsD == BnprocsD ) );

        if( !( AisR ) )
        {
            /*
            *  sub( A ) is distributed but not replicated
            */
            if( BisR )
            {
                /*
                *  If sub( A ) is not replicated, but sub( B ) is, a process row or column
                *  BprocR need to be selected. It will contain the non-replicated vector to
                *  add sub( A ) to.
                */
                if( RRorCC )
                {
                    /*
                    *  sub( A ) and sub( B ) are both row or column vectors
                    */
                    if( ( OneDgrid || OneBlock || Square ) && ( AprocD == BprocD ) )
                    {
                        /*
                        *  sub( A ) and sub( B ) start in the same process row or column AprocD=BprocD.
                        *  Enforce a purely local operation by choosing BprocR to be equal to AprocR.
                        */
                        BprocR = AprocR;
                    }
                    else
                    {
                        /*
                        *  Otherwise, communication has to occur, so choose the next process row or
                        *  column for BprocR to maximize the number of links, i.e reduce contention.
                        */
                        BprocR = MModAdd1( AprocR, AnprocsR );
                    }
                }
                else
                {
                    /*
                    *  sub( A ) and sub( B ) are distributed in orthogonal directions, what is
                    *  chosen for YprocR does not really matter. Select the process origin.
                    */
                    BprocR = AprocD;
                }
            }
            else
            {
                /*
                *  Neither sub( A ) nor sub( B ) are replicated. If I am not in process row or
                *  column AprocR and not in process row or column BprocR, then quick return.
                */
                if( ( AmyprocR != AprocR ) && ( BmyprocR != BprocR ) )
                {
                    return;
                }
            }
        }
        else
        {
            /*
            *  sub( A ) is distributed and replicated (so no quick return possible)
            */
            if( BisR )
            {
                /*
                *  sub( B ) is distributed and replicated as well
                */
                if( RRorCC )
                {
                    /*
                    *  sub( A ) and sub( B ) are both row or column vectors
                    */
                    if( ( OneDgrid || OneBlock || Square ) && ( AprocD == BprocD ) )
                    {
                        /*
                        *  sub( A ) and sub( B ) start in the same process row or column AprocD=BprocD.
                        *  Enforce a purely local operation by choosing AprocR and BprocR to be equal
                        *  to zero.
                        */
                        AprocR = BprocR = 0;
                    }
                    else
                    {
                        /*
                        *  Otherwise, communication has to occur, so select BprocR to be zero and the
                        *  next process row or column for AprocR in order to maximize the number of
                        *  used links, i.e reduce contention.
                        */
                        BprocR = 0;
                        AprocR = MModAdd1( BprocR, BnprocsR );
                    }
                }
                else
                {
                    /*
                    *  sub( A ) and sub( B ) are distributed in orthogonal directions, select the
                    *  origin processes.
                    */
                    AprocR = BprocD;
                    BprocR = AprocD;
                }
            }
            else
            {
                /*
                *  sub( B ) is distributed, but not replicated
                */
                if( RRorCC )
                {
                    /*
                    *  sub( A ) and sub( B ) are both row or column vectors
                    */
                    if( ( OneDgrid || OneBlock || Square ) && ( AprocD == BprocD ) )
                    {
                        /*
                        *  sub( A ) and sub( B ) start in the same process row or column AprocD=BprocD.
                        *  Enforce a purely local operation by choosing AprocR to be equal to BprocR.
                        */
                        AprocR = BprocR;
                        if( ( AmyprocR != AprocR ) && ( BmyprocR != BprocR ) )
                        {
                            return;
                        }
                    }
                    else
                    {
                        /*
                        *  Otherwise, communication has to occur, so choose the next process row or
                        *  column for AprocR to maximize the number of links, i.e reduce contention.
                        */
                        AprocR = MModAdd1( BprocR, BnprocsR );
                    }
                }
                else
                {
                    /*
                    *  sub( A ) and sub( B ) are distributed in orthogonal directions, what is
                    *  chosen for AprocR does not really matter. Select the origin process.
                    */
                    AprocR = BprocD;
                    if( ( OneDgrid || OneBlock || Square ) &&
                        ( AmyprocR != AprocR ) && ( BmyprocR != BprocR ) )
                    {
                        return;
                    }
                }
            }
        }

        /*
        *  Even if sub( A ) and/or sub( B ) are replicated, only two process row or
        *  column are active, namely AprocR and BprocR. If any of those operands is
        *  replicated, broadcast will occur (unless there is an easy way out).
        */
        size = TYPE->size;

        /*
        *  A purely local operation occurs iff the operands start in the same process
        *  and, if either the grid is mono-dimensional or there is a single local block
        *  to be added or if both operands are aligned.
        */
        if( ( (    RRorCC   && ( AprocD == BprocD ) &&
                         ( AisR || BisR || ( AprocR == BprocR ) ) ) ||
             ( !( RRorCC ) && ( BisR || ( AprocD == BprocR ) ) &&
                              ( AisR || ( AprocR == BprocD ) ) ) ) &&
           ( OneDgrid || OneBlock || ( RRorCC && Square ) ) )
        {
            if( ( !AisR && ( AmyprocR == AprocR ) ) ||
                ( AisR && ( BisR || BmyprocR == BprocR ) ) )
            {
                AnpD = PB_Cnumroc( AnD, 0, Ainb1D, AnbD, AmyprocD, AprocD,
                                    AnprocsD );
                BnpD = PB_Cnumroc( BnD, 0, Binb1D, BnbD, BmyprocD, BprocD,
                                    BnprocsD );
                if( ( AnpD > 0 ) && ( BnpD > 0 ) )
                {
                    /*
                    *  Select the local add routine accordingly to RRorCC
                    */
                    if( RRorCC )
                    {
                        if( Mupcase( CONJUG[0] ) != CNOCONJG )
                        {
                            // add = TYPE->Fmmcadd;
                            cout << "Conjugate multiply not implemented yet" << endl;
                            MPI_Abort(MPI_COMM_WORLD, -1);
                        }
                        else
                        {
                            mult = dmmmult;
                        }
                    }
                    else
                    {
                        if( Mupcase( CONJUG[0] ) != CNOCONJG )
                        {
                            // add = TYPE->Fmmtcadd;
                            cout << "Conjugate transpose multiply not implemented yet" << endl;
                            MPI_Abort(MPI_COMM_WORLD, -1);
                        }
                        else
                        {
                            // add = TYPE->Fmmtadd;
                            cout << "Transpose multiply not implemented yet" << endl;
                            MPI_Abort(MPI_COMM_WORLD, -1);
                        }
                    }

                    /*
                    *  Local addition
                    */
                    cout << "local addition" << endl;
                    if( AisRow )
                        mult( AnR, AnpD, ALPHA, Mptr( A, Aii, Ajj, Ald, size ), Ald,
                            Mptr( B, Bii, Bjj, Bld, size ), Bld );
                    else
                        mult( AnpD, AnR, ALPHA, Mptr( A, Aii, Ajj, Ald, size ), Ald,
                            Mptr( B, Bii, Bjj, Bld, size ), Bld );
                }
            }
            if( RRorCC && AisR && BisR )
            {
                cout << "return 4" << endl;
                return;
            }
        }
        else if( ( RRorCC && OneDgrid ) || OneBlock || Square )
        {
            /*
            *  Otherwise, it may be possible to add the distributed vectors in a single
            *  message exchange iff the grid is mono-dimensional and the operands are
            *  distributed in the same direction, or there is just one block to be exchanged
            *  or if both operands are similarly distributed in their respective direction.
            */
            if( RRorCC )
            {
                if( Mupcase( CONJUG[0] ) != CNOCONJG )
                {
                    // add = TYPE->Fmmcadd;
                    cout << "Conjugate multiply not implemented yet" << endl;
                    MPI_Abort(MPI_COMM_WORLD, -1);
                }
                else
                {
                    mult = dmmmult;
                }
            }
            else
            {
                if( Mupcase( CONJUG[0] ) != CNOCONJG )
                {
                    // add = TYPE->Fmmtcadd;
                    cout << "Conjugate transpose multiply not implemented yet" << endl;
                    MPI_Abort(MPI_COMM_WORLD, -1);
                }
                else
                {
                    // add = TYPE->Fmmtadd;
                    cout << "Transpose multiply not implemented yet" << endl;
                    MPI_Abort(MPI_COMM_WORLD, -1);
                }
            }

            if( ( AisR && BisR ) || ( AmyprocR == AprocR ) )
            {

                AnpD = PB_Cnumroc( AnD, 0, Ainb1D, AnbD, AmyprocD, AprocD,
                                    AnprocsD );
                if( AnpD > 0 )
                {
                    dst = BprocD + MModSub( AmyprocD, AprocD, AnprocsD );
                    dst = MPosMod( dst, BnprocsD );
                    if( AisRow ) { ma = AnR; na = AnpD; }
                    else         { ma = AnpD; na = AnR; }
                    if( !( AisR && BisR ) )
                    {
                        if( BisRow ) { rdst = BprocR; cdst = dst; }
                        else         { rdst = dst; cdst = BprocR; }
                    }
                    else
                    {
                        if( BisRow )
                        {
                            if( !AisRow ) { rdst = AmyprocR; }
                            else { rdst = MModAdd1( BmyprocR, BnprocsR ); }
                            cdst = dst;
                        }
                        else
                        {
                            rdst = dst;
                            if( AisRow ) { cdst = AmyprocR; }
                            else { cdst = MModAdd1( BmyprocR, BnprocsR ); }
                        }
                    }

                    if( ( myrow == rdst ) && ( mycol == cdst ) )
                    {
                        mult( ma, na, ALPHA, Mptr( A, Aii, Ajj, Ald, size ), Ald,
                            Mptr( B, Bii, Bjj, Bld, size ), Bld );
                    }
                    else
                    {
                        TYPE->Cgesd2d( ctxt, ma, na, Mptr( A, Aii, Ajj, Ald, size ),
                                    Ald, rdst, cdst );
                    }
                }
            }

            if( ( AisR && BisR ) || ( BmyprocR == BprocR ) )
            {
                BnpD = PB_Cnumroc( BnD, 0, Binb1D, BnbD, BmyprocD, BprocD,
                                    BnprocsD );
                if( BnpD > 0 )
                {
                    src = AprocD + MModSub( BmyprocD, BprocD, BnprocsD );
                    src = MPosMod( src, AnprocsD );
                    if( AisRow ) { ma = BnR; na = BnpD; }
                    else         { ma = BnpD; na = BnR; }
                    if( !( AisR && BisR ) )
                    {
                        if( AisRow ) { rsrc = AprocR; csrc = src; }
                        else         { rsrc = src; csrc = AprocR; }
                    }
                    else
                    {
                        if( AisRow )
                        {
                            if( !BisRow ) { rsrc = BmyprocR; }
                            else { rsrc = MModSub1( AmyprocR, AnprocsR ); }
                            csrc = src;
                        }
                        else
                        {
                            rsrc = src;
                            if( BisRow ) { csrc = BmyprocR; }
                            else { csrc = MModSub1( AmyprocR, AnprocsR ); }
                        }
                    }

                    if( ( myrow != rsrc ) || ( mycol != csrc ) )
                    {
                        buf = PB_Cmalloc( BnpD * BnR * size );
                        TYPE->Cgerv2d( ctxt, ma, na, buf, ma, rsrc, csrc );
                        mult( ma, na, ALPHA, buf, ma, Mptr( B, Bii, Bjj, Bld, size ), Bld );
                        if( buf ) free( buf );
                    }
                }
            }
            if( AisR && BisR )
            {
                return;
            }
        }
        else
        {

            /*
            *  General case
            */
            if( RRorCC )
            {
                if( Mupcase( CONJUG[0] ) != CNOCONJG ) tran = CCONJG;
                else                                   tran = CNOTRAN;
            }
            else
            {
                if( Mupcase( CONJUG[0] ) != CNOCONJG ) tran = CCOTRAN;
                else                                   tran = CTRAN;
            }

            if( AisRow ) { ascope = CCOLUMN; ma = AnR; }
            else         { ascope = CROW;    na = AnR; }
            bscope = ( BisRow ? CCOLUMN : CROW );
            lcmb   = PB_Clcm( AnprocsD * AnbD, BnprocsD * BnbD );
            one    = TYPE->one; zero = TYPE->zero;
            gcdPQ  = PB_Cgcd( AnprocsD, BnprocsD );
            lcmPQ  = ( AnprocsD / gcdPQ ) * BnprocsD;

            for( k = 0; k < gcdPQ; k++ )
            {
                p = 0; q = k;

                for( l = 0; l < lcmPQ; l++ )
                {
                    Aroc = MModAdd( AprocD, p, AnprocsD );
                    Broc = MModAdd( BprocD, q, BnprocsD );

                    if( ( AmyprocD == Aroc ) || ( BmyprocD == Broc ) )
                    {
                        AnpD = PB_Cnumroc( AnD, 0, Ainb1D, AnbD, Aroc, AprocD,
                                            AnprocsD );
                        BnpD = PB_Cnumroc( BnD, 0, Binb1D, BnbD, Broc, BprocD,
                                            BnprocsD );
                        PB_CVMinit( &VM, 0, AnpD, BnpD, Ainb1D, Binb1D, AnbD, BnbD,
                                    p, q, AnprocsD, BnprocsD, lcmb );
                        if( npq = PB_CVMnpq( &VM ) )
                        {
                            if( ( RRorCC && ( Aroc == Broc ) &&
                                    ( AisR || ( AprocR == BprocR ) ) ) ||
                                ( !( RRorCC ) && ( Aroc == BprocR ) &&
                                    ( AisR || ( AprocR == Broc ) ) ) )
                            {
                                if( ( BmyprocD ==  Broc ) && ( BmyprocR == BprocR ) )
                                {
                                    VMlocMult( TYPE, &VM, ROW, &ascope, PACKING, &tran,
                                            npq, AnR, ALPHA, Mptr( A, Aii, Ajj, Ald,
                                            size ), Ald, Mptr( B, Bii, Bjj, Bld,
                                            size ), Bld );
                                }
                            }
                            else
                            {
                                // send
                                if( ( AmyprocR == AprocR ) && ( AmyprocD == Aroc  ) )
                                {
                                    if( AisRow ) { na = npq; } else { ma = npq; }
                                    buf = PB_Cmalloc( ma * na * size );

                                    // pack (no addition)
                                    PB_CVMpack( TYPE, &VM, ROW, &ascope, PACKING, NOTRAN,
                                                npq, AnR, one, Mptr( A, Aii, Ajj, Ald,
                                                size ), Ald, zero, buf, ma );
                                    if( BisRow ) { rdst = BprocR; cdst = Broc; }
                                    else         { rdst = Broc; cdst = BprocR; }
                                    TYPE->Cgesd2d( ctxt, ma, na, buf, ma, rdst, cdst );
                                    if( buf ) free ( buf );
                                }
                                // recieve
                                if( ( BmyprocR == BprocR ) && ( BmyprocD == Broc ) )
                                {
                                    if( AisRow )
                                    { na = npq; rsrc = AprocR; csrc = Aroc; }
                                    else
                                    { ma = npq; rsrc = Aroc; csrc = AprocR; }
                                    buf = PB_Cmalloc( ma * na * size );

                                    // receive from Cgesd2d
                                    TYPE->Cgerv2d( ctxt, ma, na, buf, ma, rsrc, csrc );

                                    // recieve, then multiply
                                    VMpackMult( TYPE, &VM, COLUMN, &bscope, UNPACKING,
                                                &tran, npq, AnR, ALPHA, Mptr( B, Bii, Bjj,
                                                Bld, size ), Bld, buf, ma );
                                    if( buf ) free ( buf );
                                }
                            }
                        }
                    }
                    p = MModAdd1( p, AnprocsD );
                    q = MModAdd1( q, BnprocsD );
                }
                if( AisR ) AprocR = MModAdd1( AprocR, AnprocsR );
            }
        }

        if( BisR )
        {
            /*
            *  Replicate sub( B )
            */
            BnpD = PB_Cnumroc( BnD, BiD, BinbD, BnbD, BmyprocD, BsrcD, BnprocsD );
            if( BnpD > 0 )
            {
                if( BisRow )
                {
                    bscope = CCOLUMN;  mb   = BnR;      nb = BnpD;
                    rsrc   = BprocR;   csrc = BmyprocD;
                }
                else
                {
                    bscope = CROW;     mb   = BnpD;     nb = BnR;
                    rsrc   = BmyprocD; csrc = BprocR;
                }
                top = PB_Ctop( &ctxt, BCAST, &bscope, TOP_GET );
                if( BmyprocR == BprocR )
                {
                    TYPE->Cgebs2d( ctxt, &bscope, top, mb, nb, Mptr( B, Bii, Bjj,
                                Bld, size ), Bld );
                }
                else
                {
                    TYPE->Cgebr2d( ctxt, &bscope, top, mb, nb, Mptr( B, Bii, Bjj,
                                Bld, size ), Bld, rsrc, csrc );
                }
            }
        }
    }
    else
    {
        cout << "daxy() cannot handle mismatched distribution" << endl;
        cout << "You'll have to manually implement PB_CpaxpbyND/DN/NN" << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

}

int VMlocMult(PBTYP_T * TYPE, PB_VM_T * VM, char * VROCS, char * ROCS,
    char * UNPA, char * TRANS, int MN, int K, char * ALPHA,
    char * A, int LDA,
    char * B, int LDB)
{
    /*
    *  .. Local Scalars ..
    */
    int            GoEast, GoSouth, ilow, imbloc, inbloc, inca, incb, iupp, kb,
                   lcmt, lcmt00, lmbloc, lnbloc, low, mb, mblkd, mblks, mbloc,
                   * m, * n, nb, nblkd, nblks, nbloc, notran, npcol, npq=0,
                   nprow, pmb, qnb, rows, size, tmp1, tmp2, upp;
    MMMULT_T       mult;
    char           * aptrd, * bptrd;

    MMMULT_T dmmmult = dmmmult_;
    MMMULT_T dmmtlum = dmmtlum_;

    mblks = VM->mblks; nblks = VM->nblks;

    /*
    *  Quick return if I don't own any blocks.
    */
    if( ( mblks == 0 ) || ( nblks == 0 ) ) return( 0 );

    /*
    *  Retrieve the contents of VM structure fields
    */
    lcmt00 = VM->lcmt00;
    imbloc = VM->imbloc; mb    = VM->mb; lmbloc = VM->lmbloc; upp = VM->upp;
    iupp   = VM->iupp;   nprow = VM->nprow;
    inbloc = VM->inbloc; nb    = VM->nb; lnbloc = VM->lnbloc; low = VM->low;
    ilow   = VM->ilow;   npcol = VM->npcol;

    if( Mupcase( UNPA[0] ) == CPACKING )
    {
        /*
        *  B is the distributed target, A is the distributed source
        */
        if( Mupcase( TRANS[0] ) == CNOTRAN )
        {
            /*
            *  Add A to B
            */
            notran = 1;
            mult = dmmmult;
        }
        else if( Mupcase( TRANS[0] ) == CCONJG )
        {
            /*
            *  Add the conjugate of A to B
            */
            // notran = 1;
            // add    = TYPE->Fmmcadd;
            cout << "Conjugate multiply is not implemented yet" << endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        else if( Mupcase( TRANS[0] ) == CTRAN )
        {
            /*
            *  Add the tranpose of A to B
            */
            // notran = 0;
            // add    = TYPE->Fmmtadd;
            cout << "Transpose multiply is not implemented yet" << endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        else
        {
            /*
            *  Add the conjugate tranpose of A to B
            */
            // notran = 0;
            // add    = TYPE->Fmmtcadd;
            cout << "Conjugate transpose multiply is not implemented yet" << endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }
    else
    {
        /*
        *  A is the distributed target, B is the distributed source
        */
        if( Mupcase( TRANS[0] ) == CNOTRAN )
        {
            /*
            *  Add B to A
            */
            notran = 1;
            mult = dmmtlum;
        }
        else if( Mupcase( TRANS[0] ) == CCONJG )
        {
            /*
            *  Add the conjugate of B to A
            */
            // notran = 1;
            // add    = TYPE->Fmmddac;
            cout << "Conjugate multiply is not implemented yet" << endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        else if( Mupcase( TRANS[0] ) == CTRAN )
        {
            /*
            *  Add the tranpose of B to A
            */
            // notran = 0;
            // add    = TYPE->Fmmddat;
            cout << "Transpose multiply is not implemented yet" << endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        else
        {
            /*
            *  Add the conjugate tranpose of B to A
            */
            // notran = 0;
            // add    = TYPE->Fmmddact;
            cout << "Conjugate tranpose multiply is not implemented yet" << endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }

    size = TYPE->size;
    rows = ( Mupcase( ROCS[0] ) == CROW );

    if( Mupcase( VROCS[0] ) == CROW )
    {
        /*
        *  (un)packing using rows of virtual matrix
        */
        if( rows )
        {
            /*
            *  (un)packing rows of mn by k array A.
            */
            inca = size;
            incb = ( notran ? size : LDB * size );
            m    = &tmp2;
            n    = &K;
        }
        else
        {
            /*
            *  (un)packing columns of k by mn array A
            */
            inca = LDA * size;
            incb = ( notran ? LDB * size : size );
            m    = &K;
            n    = &tmp2;
        }
        kb  = MN;

        /*
        *  From the (un)packing point of view the only valuable shortcut is when the
        *  virtual grid and the blocks are square, and the offset is zero or the grid
        *  is 1x1.
        */
        if( ( ( lcmt00 == 0 ) && ( VM->imb1 == VM->inb1 ) && ( mb == nb ) &&
                ( nprow == npcol ) ) || ( ( nprow == 1 ) && ( npcol == 1 ) ) )
        {
            if(  VM->prow == VM->pcol )
            {
                npq = ( ( mblks <  2 ) ? imbloc :
                        imbloc + ( mblks - 2 ) * mb + lmbloc );
                npq = min( npq, kb );
                if( rows ) mult( npq, K, ALPHA, A, LDA, B, LDB );
                else       mult( K, npq, ALPHA, A, LDA, B, LDB );
            }
            return( npq );
        }
        pmb = nprow * mb;
        qnb = npcol * nb;

        /*
        *  Handle separately the first row and/or column of the LCM table. Update the
        *  LCM value of the curent block lcmt00, as well as the number of rows and
        *  columns mblks and nblks remaining in the LCM table.
        */
        GoSouth = ( lcmt00 > iupp );
        GoEast  = ( lcmt00 < ilow );

        /*
        *  Go through the table looking for blocks owning diagonal entries.
        */
        if( !( GoSouth ) && !( GoEast ) )
        {

            /*
            *  The upper left block owns diagonal entries lcmt00 >= ilow && lcmt00 <= iupp
            */
            if( lcmt00 >= 0 )
            {
                tmp1 = imbloc - lcmt00; tmp1 = max( 0, tmp1 );
                tmp2 = min( tmp1, inbloc ); npq += ( tmp2 = min( tmp2, kb ) );
                mult( *m, *n, ALPHA, A+lcmt00*inca, LDA, B, LDB );
            }
            else
            {
                tmp1 = inbloc + lcmt00; tmp1 = max( 0, tmp1 );
                tmp2 = min( tmp1, imbloc ); npq += ( tmp2 = min( tmp2, kb ) );
                mult( *m, *n, ALPHA, A, LDA, B-lcmt00*incb, LDB );
            }
            if( ( kb -= tmp2 ) == 0 ) return( npq );
            /*
            *  Decide whether one should go south or east in the table: Go east if
            *  the block below the current one only owns lower entries. If this block,
            *  however, owns diagonals, then go south.
            */
            GoSouth = !( GoEast = ( ( lcmt00 - ( iupp - upp + pmb ) ) < ilow ) );
        }

        if( GoSouth )
        {
            /*
            *  Go one step south in the LCM table. Adjust the current LCM value as well as
            *  the pointer to A. The pointer to B remains unchanged.
            */
            lcmt00 -= iupp - upp + pmb; mblks--; A += imbloc * inca;

            /*
            *  While there are blocks remaining that own upper entries, keep going south.
            *  Adjust the current LCM value as well as the pointer to A accordingly.
            */
            while( mblks && ( lcmt00 > upp ) )
            { lcmt00 -= pmb; mblks--; A += mb * inca; }

            /*
            *  Return if no more row in the LCM table.
            */
            if( mblks <= 0 ) return( npq );

            /*
            *  lcmt00 <= upp. The current block owns either diagonals or lower entries.
            *  Save the current position in the LCM table. After this column has been
            *  completely taken care of, re-start from this row and the next column of
            *  the LCM table.
            */
            lcmt = lcmt00; mblkd = mblks; aptrd = A;

            while( mblkd && ( lcmt >= ilow ) )
            {
                /*
                *  A block owning diagonals lcmt00 >= ilow && lcmt00 <= upp has been found.
                */
                mbloc = ( ( mblkd == 1 ) ? lmbloc : mb );
                if( lcmt >= 0 )
                {
                    tmp1 = mbloc - lcmt; tmp1 = max( 0, tmp1 );
                    tmp2 = min( tmp1, inbloc ); npq += ( tmp2 = min( tmp2, kb ) );
                    mult( *m, *n, ALPHA, aptrd+lcmt*inca, LDA, B, LDB );
                }
                else
                {
                    tmp1 = inbloc + lcmt; tmp1 = max( 0, tmp1 );
                    tmp2 = min( tmp1, mbloc ); npq += ( tmp2 = min( tmp2, kb ) );
                    mult( *m, *n, ALPHA, aptrd, LDA, B-lcmt*incb, LDB );
                }
                if( ( kb -= tmp2 ) == 0 ) return( npq );
                /*
                *  Keep going south until there are no more blocks owning diagonals
                */
                lcmt -= pmb; mblkd--; aptrd += mbloc * inca;
            }

            /*
            *  I am done with the first column of the LCM table. Go to the next column.
            */
            lcmt00 += low - ilow + qnb; nblks--; B += inbloc * incb;
        }
        else if( GoEast )
        {
            /*
            *  Go one step east in the LCM table. Adjust the current LCM value as
            *  well as the pointer to B. The pointer to A remains unchanged.
            */
            lcmt00 += low - ilow + qnb; nblks--; B += inbloc * incb;

            /*
            *  While there are blocks remaining that own lower entries, keep going east
            *  in the LCM table. Adjust the current LCM value as well as the pointer to
            *  B accordingly.
            */
            while( nblks && ( lcmt00 < low ) )
            { lcmt00 += qnb; nblks--; B += nb * incb; }

            /*
            *  Return if no more column in the LCM table.
            */
            if( nblks <= 0 ) return( npq );

            /*
            *  lcmt00 >= low. The current block owns either diagonals or upper entries. Save
            *  the current position in the LCM table. After this row has been completely
            *  taken care of, re-start from this column and the next row of the LCM table.
            */
            lcmt = lcmt00; nblkd = nblks; bptrd = B;

            while( nblkd && ( lcmt <= iupp ) )
            {
                /*
                *  A block owning diagonals lcmt00 >= low && lcmt00 <= iupp has been found.
                */
                nbloc = ( ( nblkd == 1 ) ? lnbloc : nb );
                if( lcmt >= 0 )
                {
                    tmp1 = imbloc - lcmt; tmp1 = max( 0, tmp1 );
                    tmp2 = min( tmp1, nbloc ); npq += ( tmp2 = min( tmp2, kb ) );
                    mult( *m, *n, ALPHA, A+lcmt*inca, LDA, bptrd, LDB );
                }
                else
                {
                    tmp1 = nbloc + lcmt; tmp1 = max( 0, tmp1 );
                    tmp2 = min( tmp1, imbloc ); npq += ( tmp2 = min( tmp2, kb ) );
                    mult( *m, *n, ALPHA, A, LDA, bptrd-lcmt*incb, LDB );
                }
                if( ( kb -= tmp2 ) == 0 ) return( npq );

                /*
                *  Keep going east until there are no more blocks owning diagonals.
                */
                lcmt += qnb; nblkd--; bptrd += nbloc * incb;
            }

            /*
            *  I am done with the first row of the LCM table. Go to the next row.
            */
            lcmt00 -= iupp - upp + pmb; mblks--; A += imbloc * inca;
        }
        /*
        *  Loop over the remaining columns of the LCM table.
        */
        do
        {
            /*
            *  If the current block does not have diagonal elements, find the closest one in
            *  the LCM table having some.
            */
            if( ( lcmt00 < low ) || ( lcmt00 > upp ) )
            {
                while( mblks && nblks )
                {
                    while( mblks && ( lcmt00 > upp ) )
                    { lcmt00 -= pmb; mblks--; A += mb * inca; }
                    if( lcmt00 >= low ) break;
                    while( nblks && ( lcmt00 < low ) )
                    { lcmt00 += qnb; nblks--; B += nb * incb; }
                    if( lcmt00 <= upp ) break;
                }
            }
            if( !( mblks ) || !( nblks ) ) return( npq );

            /*
            *  The current block owns diagonals. Save the current position in the LCM table.
            *  After this column has been completely taken care of, re-start from this row
            *  and the next column in the LCM table.
            */
            nbloc = ( ( nblks == 1 ) ? lnbloc : nb );
            lcmt = lcmt00; mblkd = mblks; aptrd = A;

            while( mblkd && ( lcmt >= low ) )
            {
                /*
                *  A block owning diagonals lcmt00 >= low && lcmt00 <= upp has been found.
                */
                mbloc = ( ( mblkd == 1 ) ? lmbloc : mb );
                if( lcmt >= 0 )
                {
                    tmp1 = mbloc - lcmt; tmp1 = max( 0, tmp1 );
                    tmp2 = min( tmp1, nbloc ); npq += ( tmp2 = min( tmp2, kb ) );
                    mult( *m, *n, ALPHA, aptrd+lcmt*inca, LDA, B, LDB );
                }
                else
                {
                    tmp1 = nbloc + lcmt; tmp1 = max( 0, tmp1 );
                    tmp2 = min( tmp1, mbloc ); npq += ( tmp2 = min( tmp2, kb ) );
                    mult( *m, *n, ALPHA, aptrd, LDA, B-lcmt*incb, LDB );
                }
                if( ( kb -= tmp2 ) == 0 ) return( npq );

                /*
                *  Keep going south until there are no more blocks owning diagonals
                */
                lcmt -= pmb; mblkd--; aptrd += mbloc * inca;
            }
            /*
            *  I am done with this column of the LCM table. Go to the next column ...
            */
            lcmt00 += qnb; nblks--; B += nbloc * incb;

            /*
            *  ... until there are no more columns.
            */
        } while( nblks > 0 );

        /*
        *  Return the number of diagonals found.
        */
        return( npq );
    }
    else
    {
        /*
        *  (un)packing using columns of virtual matrix
        */
        if( rows )
        {
            /*
            *  (un)packing rows of mn by k array A
            */
            inca = size;
            incb = ( notran ? size : LDB * size );
            m    = &tmp2;
            n    = &K;
        }
        else
        {
            /*
            *  (un)packing columns of k by mn array A
            */
            inca = LDA * size;
            incb = ( notran ? LDB * size : size );
            m    = &K;
            n    = &tmp2;
        }
        kb  = MN;

        /*
        *  From the (un)packing point of view the only valuable shortcut is when the
        *  virtual grid and the blocks are square, and the offset is zero or the grid
        *  is 1x1.
        */
        if( ( ( lcmt00 == 0 ) && ( VM->imb1 == VM->inb1 ) && ( mb == nb ) &&
                ( nprow == npcol ) ) || ( ( nprow == 1 ) && ( npcol == 1 ) ) )
        {
            if(  VM->prow == VM->pcol )
            {
                npq = ( ( nblks <  2 ) ? inbloc :
                        inbloc + ( nblks - 2 ) * nb + lnbloc );
                npq = min( npq, kb );
                if( rows ) mult( npq, K, ALPHA, A, LDA, B, LDB );
                else       mult( K, npq, ALPHA, A, LDA, B, LDB );
            }
            return( npq );
        }
        pmb = nprow * mb;
        qnb = npcol * nb;

        /*
        *  Handle separately the first row and/or column of the LCM table. Update the
        *  LCM value of the curent block lcmt00, as well as the number of rows and
        *  columns mblks and nblks remaining in the LCM table.
        */
        GoSouth = ( lcmt00 > iupp );
        GoEast  = ( lcmt00 < ilow );

        if( !( GoSouth ) && !( GoEast ) )
        {
            /*
            *  The upper left block owns diagonal entries lcmt00 >= ilow && lcmt00 <= iupp
            */
            if( lcmt00 >= 0 )
            {
                tmp1 = imbloc - lcmt00; tmp1 = max( 0, tmp1 );
                tmp2 = min( tmp1, inbloc ); npq += ( tmp2 = min( tmp2, kb ) );
                mult( *m, *n, ALPHA, A, LDA, B+lcmt00*incb, LDB );
            }
            else
            {
                tmp1 = inbloc + lcmt00; tmp1 = max( 0, tmp1 );
                tmp2 = min( tmp1, imbloc ); npq += ( tmp2 = min( tmp2, kb ) );
                mult( *m, *n, ALPHA, A-lcmt00*inca, LDA, B, LDB );
            }
            if( ( kb -= tmp2 ) == 0 ) return( npq );

            /*
            *  Decide whether one should go south or east in the table: Go east if
            *  the block below the current one only owns lower entries. If this block,
            *  however, owns diagonals, then go south.
            */
            GoSouth = !( GoEast = ( ( lcmt00 - ( iupp - upp + pmb ) ) < ilow ) );
        }

        if( GoSouth )
        {
            /*
            *  Go one step south in the LCM table. Adjust the current LCM value as well as
            *  the pointer to B. The pointer to A remains unchanged.
            */
            lcmt00 -= iupp - upp + pmb; mblks--; B += imbloc * incb;

            /*
            *  While there are blocks remaining that own upper entries, keep going south.
            *  Adjust the current LCM value as well as the pointer to B accordingly.
            */
            while( mblks && ( lcmt00 > upp ) )
            { lcmt00 -= pmb; mblks--; B += mb * incb; }

            /*
            *  Return if no more row in the LCM table.
            */
            if( mblks <= 0 ) return( npq );

            /*
            *  lcmt00 <= upp. The current block owns either diagonals or lower entries.
            *  Save the current position in the LCM table. After this column has been
            *  completely taken care of, re-start from this row and the next column of
            *  the LCM table.
            */
            lcmt  = lcmt00; mblkd = mblks; bptrd = B;

            while( mblkd && ( lcmt >= ilow ) )
            {
                /*
                *  A block owning diagonals lcmt00 >= ilow && lcmt00 <= upp has been found.
                */
                mbloc = ( ( mblkd == 1 ) ? lmbloc : mb );
                if( lcmt >= 0 )
                {
                    tmp1 = mbloc - lcmt; tmp1 = max( 0, tmp1 );
                    tmp2 = min( tmp1, inbloc ); npq += ( tmp2 = min( tmp2, kb ) );
                    mult( *m, *n, ALPHA, A, LDA, bptrd+lcmt*incb, LDB );
                }
                else
                {
                    tmp1 = inbloc + lcmt; tmp1 = max( 0, tmp1 );
                    tmp2 = min( tmp1, mbloc ); npq += ( tmp2 = min( tmp2, kb ) );
                    mult( *m, *n, ALPHA, A-lcmt*inca, LDA, bptrd, LDB );
                }
                if( ( kb -= tmp2 ) == 0 ) return( npq );

                /*
                *  Keep going south until there are no more blocks owning diagonals
                */
                lcmt -= pmb; mblkd--; bptrd += mbloc * incb;
            }

            /*
            *  I am done with the first column of the LCM table. Go to the next column.
            */
            lcmt00 += low - ilow + qnb; nblks--; A += inbloc * inca;
        }
        else if( GoEast )
        {
            /*
            *  Go one step east in the LCM table. Adjust the current LCM value as
            *  well as the pointer to A. The pointer to B remains unchanged.
            */
            lcmt00 += low - ilow + qnb; nblks--; A += inbloc * inca;

            /*
            *  While there are blocks remaining that own lower entries, keep going east
            *  in the LCM table. Adjust the current LCM value as well as the pointer to
            *  A accordingly.
            */
            while( nblks && ( lcmt00 < low ) )
            { lcmt00 += qnb; nblks--; A += nb * inca; }

            /*
            *  Return if no more column in the LCM table.
            */
            if( nblks <= 0 ) return( npq );

            /*
            *  lcmt00 >= low. The current block owns either diagonals or upper entries. Save
            *  the current position in the LCM table. After this row has been completely
            *  taken care of, re-start from this column and the next row of the LCM table.
            */
            lcmt  = lcmt00; nblkd = nblks; aptrd = A;

            while( nblkd && ( lcmt <= iupp ) )
            {
                /*
                *  A block owning diagonals lcmt00 >= low && lcmt00 <= iupp has been found.
                */
                nbloc = ( ( nblkd == 1 ) ? lnbloc : nb );
                if( lcmt >= 0 )
                {
                    tmp1 = imbloc - lcmt; tmp1 = max( 0, tmp1 );
                    tmp2 = min( tmp1, nbloc ); npq += ( tmp2 = min( tmp2, kb ) );
                    mult( *m, *n, ALPHA, aptrd, LDA, B+lcmt*incb, LDB );
                }
                else
                {
                    tmp1 = nbloc + lcmt; tmp1 = max( 0, tmp1 );
                    tmp2 = min( tmp1, imbloc ); npq += ( tmp2 = min( tmp2, kb ) );
                    mult( *m, *n, ALPHA, aptrd-lcmt*inca, LDA, B, LDB );
                }
                if( ( kb -= tmp2 ) == 0 ) return( npq );

                /*
                *  Keep going east until there are no more blocks owning diagonals.
                */
                lcmt += qnb; nblkd--; aptrd += nbloc * inca;
            }

            /*
            *  I am done with the first row of the LCM table. Go to the next row.
            */
            lcmt00 -= iupp - upp + pmb; mblks--; B += imbloc * incb;
        }

        /*
        *  Loop over the remaining columns of the LCM table.
        */
        do
        {
            /*
            *  If the current block does not have diagonal elements, find the closest one in
            *  the LCM table having some.
            */
            if( ( lcmt00 < low ) || ( lcmt00 > upp ) )
            {
                while( mblks && nblks )
                {
                    while( mblks && ( lcmt00 > upp ) )
                    { lcmt00 -= pmb; mblks--; B += mb * incb; }
                    if( lcmt00 >= low ) break;
                    while( nblks && ( lcmt00 < low ) )
                    { lcmt00 += qnb; nblks--; A += nb * inca; }
                    if( lcmt00 <= upp ) break;
                }
            }
            if( !( mblks ) || !( nblks ) ) return( npq );

            /*
            *  The current block owns diagonals. Save the current position in the LCM table.
            *  After this column has been completely taken care of, re-start from this row
            *  and the next column in the LCM table.
            */
            nbloc = ( ( nblks == 1 ) ? lnbloc : nb );
            lcmt = lcmt00; mblkd = mblks; bptrd = B;

            while( mblkd && ( lcmt >= low ) )
            {
                /*
                *  A block owning diagonals lcmt00 >= low && lcmt00 <= upp has been found.
                */
                mbloc = ( ( mblkd == 1 ) ? lmbloc : mb );
                if( lcmt >= 0 )
                {
                    tmp1 = mbloc - lcmt; tmp1 = max( 0, tmp1 );
                    tmp2 = min( tmp1, nbloc ); npq += ( tmp2 = min( tmp2, kb ) );
                    mult( *m, *n, ALPHA, A, LDA, bptrd+lcmt*incb, LDB );
                }
                else
                {
                    tmp1 = nbloc + lcmt; tmp1 = max( 0, tmp1 );
                    tmp2 = min( tmp1, mbloc ); npq += ( tmp2 = min( tmp2, kb ) );
                    mult( *m, *n, ALPHA, A-lcmt*inca, LDA, bptrd, LDB );
                }
                if( ( kb -= tmp2 ) == 0 ) return( npq );

                /*
                *  Keep going south until there are no more blocks owning diagonals
                */
                lcmt -= pmb; mblkd--; bptrd += mbloc * incb;
            }

            /*
            *  I am done with this column of the LCM table. Go to the next column ...
            */
            lcmt00 += qnb; nblks--; A += nbloc * inca;

            /*
            *  ... until there are no more columns.
            */
        } while( nblks > 0 );
        /*
        *  Return the number of diagonals found.
        */
        return( npq );
    }
}

int VMpackMult(PBTYP_T* TYPE, PB_VM_T* VM, char* VROCS, char* ROCS,
    char* UNPA, char* TRANS, int MN, int K, char* ALPHA,
    char* A, int LDA,
    char* B, int LDB)
{

    /*
    *  .. Local Scalars ..
    */
    int            GoEast, GoSouth, ilow, imbloc, inbloc, inca, incb, iupp, kb,
                   lcmt, lcmt00, lmbloc, lnbloc, low, mb, mblkd, mblks, mbloc,
                   * m, * n, nb, nblkd, nblks, nbloc, notran, npcol, npq=0,
                   nprow, pmb, qnb, rows, size, tmp1, tmp2, upp;
    char           * aptrd;
    MMMULT_T       mult;


    mblks = VM->mblks; nblks = VM->nblks;

    /*
    *  Quick return if I don't own any blocks.
    */
    if( ( mblks == 0 ) || ( nblks == 0 ) ) return( 0 );

    /*
    *  Retrieve the contents of VM structure fields
    */
    lcmt00 = VM->lcmt00;
    imbloc = VM->imbloc; mb    = VM->mb; lmbloc = VM->lmbloc; upp = VM->upp;
    iupp   = VM->iupp;   nprow = VM->nprow;
    inbloc = VM->inbloc; nb    = VM->nb; lnbloc = VM->lnbloc; low = VM->low;
    ilow   = VM->ilow;   npcol = VM->npcol;

    MMMULT_T dmmmult = dmmmult_;
    MMMULT_T dmmtlum = dmmtlum_;

    if( Mupcase( UNPA[0] ) == CPACKING )
    {
        /*
        *  B is the target packed buffer, A is the distributed source
        */
        if( Mupcase( TRANS[0] ) == CNOTRAN )
        {
            /*
            *  Add A to B
            */
            notran = 1;
            mult = dmmmult;
        }
        else if( Mupcase( TRANS[0] ) == CCONJG )
        {
            /*
            *  Add the conjugate of A to B
            */
            // notran = 1; add = TYPE->Fmmcadd;
            cout << "Conjugate multiply not working yet" << endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        else if( Mupcase( TRANS[0] ) == CTRAN )
        {
            /*
            *  Add the tranpose of A to B
            */
            // notran = 0; add = TYPE->Fmmtadd;
            cout << "Transpose multiply not working yet" << endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        else
        {
            /*
            *  Add the conjugate tranpose of A to B
            */
            // notran = 0; add = TYPE->Fmmtcadd;
            cout << "Conjugate transpose multiply not working yet" << endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }
    else
    {
        /*
        *  B is the source packed buffer, A is the distributed target
        */
        if( Mupcase( TRANS[0] ) == CNOTRAN )
        {
            /*
            *  Add B to A
            */
            notran = 1;
            mult = dmmtlum;
        }
        else if( Mupcase( TRANS[0] ) == CCONJG )
        {
            /*
            *  Add the conjugate of B to A
            */
            // notran = 1; add = TYPE->Fmmddac;
            cout << "Conjugate multiply not working yet" << endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        else if( Mupcase( TRANS[0] ) == CTRAN )
        {
            /*
            *  Add the tranpose of B to A
            */
            // notran = 0; add = TYPE->Fmmddat;
            cout << "Transpose multiply not working yet" << endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        else
        {
            /*
            *  Add the conjugate tranpose of B to A
            */
            // notran = 0; add = TYPE->Fmmddact;
            cout << "Conjugate transpose multiply not working yet" << endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }

    size = TYPE->size;
    rows = ( Mupcase( ROCS[0] ) == CROW );

    if( Mupcase( VROCS[0] ) == CROW )
    {
        /*
        *  (un)packing using rows of virtual matrix
        */
        if( rows )
        {
            /*
            *  (un)packing rows of mn by k array A.
            */
            inca = size;
            incb = ( notran ? size : LDB * size );
            m    = &tmp2;
            n    = &K;
        }
        else
        {
            /*
            *  (un)packing columns of k by mn array A
            */
            inca = LDA * size;
            incb = ( notran ? LDB * size : size );
            m    = &K;
            n    = &tmp2;
        }
        kb  = MN;

        /*
        *  From the (un)packing point of view the only valuable shortcut is when the
        *  virtual grid and the blocks are square, and the offset is zero or the grid
        *  is 1x1.
        */
        if( ( ( lcmt00 == 0 ) && ( VM->imb1 == VM->inb1 ) && ( mb == nb ) &&
                ( nprow == npcol ) ) || ( ( nprow == 1 ) && ( npcol == 1 ) ) )
        {
            if( VM->prow == VM->pcol )
            {
                npq = ( ( mblks <  2 ) ? imbloc :
                        imbloc + ( mblks - 2 ) * mb + lmbloc );
                npq = min( npq, kb );
                if( rows ) mult( npq, K, ALPHA, A, LDA, B, LDB );
                else       mult( K, npq, ALPHA, A, LDA, B, LDB );
            }
            return( npq );
        }
        pmb = nprow * mb;
        qnb = npcol * nb;

        /*
        *  Handle separately the first row and/or column of the LCM table. Update the
        *  LCM value of the curent block lcmt00, as well as the number of rows and
        *  columns mblks and nblks remaining in the LCM table.
        */
        GoSouth = ( lcmt00 > iupp );
        GoEast  = ( lcmt00 < ilow );

        if( !( GoSouth ) && !( GoEast ) )
        {
            /*
            *  The upper left block owns diagonal entries lcmt00 >= ilow && lcmt00 <= iupp
            */
            if( lcmt00 >= 0 )
            {
                tmp1 = imbloc - lcmt00; tmp1 = max( 0, tmp1 );
                tmp2 = min( tmp1, inbloc ); npq += ( tmp2 = min( tmp2, kb ) );
                mult( *m, *n, ALPHA, A+lcmt00*inca, LDA, B, LDB );
            }
            else
            {
                tmp1 = inbloc + lcmt00; tmp1 = max( 0, tmp1 );
                tmp2 = min( tmp1, imbloc ); npq += ( tmp2 = min( tmp2, kb ) );
                mult( *m, *n, ALPHA, A, LDA, B, LDB );
            }
            if( ( kb -= tmp2 ) == 0 ) return( npq );
            B += tmp2 * incb;

            /*
            *  Decide whether one should go south or east in the table: Go east if
            *  the block below the current one only owns lower entries. If this block,
            *  however, owns diagonals, then go south.
            */
            GoSouth = !( GoEast = ( ( lcmt00 - ( iupp - upp + pmb ) ) < ilow ) );
        }

        if( GoSouth )
        {
            /*
            *  Go one step south in the LCM table. Adjust the current LCM value as well as
            *  the pointer to A. The pointer to B remains unchanged.
            */
            lcmt00 -= iupp - upp + pmb; mblks--; A += imbloc * inca;

            /*
            *  While there are blocks remaining that own upper entries, keep going south.
            *  Adjust the current LCM value as well as the pointer to A accordingly.
            */
            while( mblks && ( lcmt00 > upp ) )
            { lcmt00 -= pmb; mblks--; A += mb * inca; }

            /*
            *  Return if no more row in the LCM table.
            */
            if( mblks <= 0 ) return( npq );

            /*
            *  lcmt00 <= upp. The current block owns either diagonals or lower entries.
            *  Save the current position in the LCM table. After this column has been
            *  completely taken care of, re-start from this row and the next column of
            *  the LCM table.
            */
            lcmt = lcmt00; mblkd = mblks; aptrd = A;

            while( mblkd && ( lcmt >= ilow ) )
            {
                /*
                *  A block owning diagonals lcmt00 >= ilow && lcmt00 <= upp has been found.
                */
                mbloc = ( ( mblkd == 1 ) ? lmbloc : mb );
                if( lcmt >= 0 )
                {
                    tmp1 = mbloc - lcmt; tmp1 = max( 0, tmp1 );
                    tmp2 = min( tmp1, inbloc ); npq += ( tmp2 = min( tmp2, kb ) );
                    mult( *m, *n, ALPHA, aptrd+lcmt*inca, LDA, B, LDB );
                }
                else
                {
                    tmp1 = inbloc + lcmt; tmp1 = max( 0, tmp1 );
                    tmp2 = min( tmp1, mbloc ); npq += ( tmp2 = min( tmp2, kb ) );
                    mult( *m, *n, ALPHA, aptrd, LDA, B, LDB );
                }
                if( ( kb -= tmp2 ) == 0 ) return( npq );

                /*
                *  Keep going south until there are no more blocks owning diagonals
                */
                lcmt -= pmb; mblkd--; aptrd += mbloc * inca; B += tmp2 * incb;
            }

            /*
            *  I am done with the first column of the LCM table. Go to the next column.
            */
            lcmt00 += low - ilow + qnb; nblks--;
        }
        else if( GoEast )
        {
            /*
            *  Go one step east in the LCM table. Adjust the current LCM value as
            *  well as the pointer to B. The pointer to A remains unchanged.
            */
            lcmt00 += low - ilow + qnb; nblks--;

            /*
            *  While there are blocks remaining that own lower entries, keep going east
            *  in the LCM table. Adjust the current LCM value as well as the pointer to
            *  B accordingly.
            */
            while( nblks && ( lcmt00 < low ) ) { lcmt00 += qnb; nblks--; }

            /*
            *  Return if no more column in the LCM table.
            */
            if( nblks <= 0 ) return( npq );

            /*
            *  lcmt00 >= low. The current block owns either diagonals or upper entries. Save
            *  the current position in the LCM table. After this row has been completely
            *  taken care of, re-start from this column and the next row of the LCM table.
            */
            lcmt = lcmt00; nblkd = nblks;

            while( nblkd && ( lcmt <= iupp ) )
            {
                /*
                *  A block owning diagonals lcmt00 >= low && lcmt00 <= iupp has been found.
                */
                nbloc = ( ( nblkd == 1 ) ? lnbloc : nb );
                if( lcmt >= 0 )
                {
                    tmp1 = imbloc - lcmt; tmp1 = max( 0, tmp1 );
                    tmp2 = min( tmp1, nbloc ); npq += ( tmp2 = min( tmp2, kb ) );
                    mult( *m, *n, ALPHA, A+lcmt*inca, LDA, B, LDB );
                }
                else
                {
                    tmp1 = nbloc + lcmt; tmp1 = max( 0, tmp1 );
                    tmp2 = min( tmp1, imbloc ); npq += ( tmp2 = min( tmp2, kb ) );
                    mult( *m, *n, ALPHA, A, LDA, B, LDB );
                }
                if( ( kb  -= tmp2 ) == 0 ) return( npq );

                /*
                *  Keep going east until there are no more blocks owning diagonals.
                */
                lcmt += qnb; nblkd--; B += tmp2 * incb;
            }

            /*
            *  I am done with the first row of the LCM table. Go to the next row.
            */
            lcmt00 -= iupp - upp + pmb; mblks--; A += imbloc * inca;
        }

        /*
        *  Loop over the remaining columns of the LCM table.
        */
        do
        {
            /*
            *  If the current block does not have diagonal elements, find the closest one in
            *  the LCM table having some.
            */
            if( ( lcmt00 < low ) || ( lcmt00 > upp ) )
            {
                while( mblks && nblks )
                {
                    while( mblks && ( lcmt00 > upp ) )
                    { lcmt00 -= pmb; mblks--; A += mb*inca; }
                    if( lcmt00 >= low ) break;
                    while( nblks && ( lcmt00 < low ) )
                    { lcmt00 += qnb; nblks--; }
                    if( lcmt00 <= upp ) break;
                }
            }
            if( !( mblks ) || !( nblks ) ) return( npq );

            /*
            *  The current block owns diagonals. Save the current position in the LCM table.
            *  After this column has been completely taken care of, re-start from this row
            *  and the next column in the LCM table.
            */
            nbloc = ( ( nblks == 1 ) ? lnbloc : nb );
            lcmt = lcmt00; mblkd = mblks; aptrd = A;

            while( mblkd && ( lcmt >= low ) )
            {
                /*
                *  A block owning diagonals lcmt00 >= low && lcmt00 <= upp has been found.
                */
                mbloc = ( ( mblkd == 1 ) ? lmbloc : mb );
                if( lcmt >= 0 )
                {
                    tmp1 = mbloc - lcmt; tmp1 = max( 0, tmp1 );
                    tmp2 = min( tmp1, nbloc ); npq += ( tmp2 = min( tmp2, kb ) );
                    mult( *m, *n, ALPHA, aptrd+lcmt*inca, LDA, B, LDB );
                }
                else
                {
                    tmp1 = nbloc + lcmt; tmp1 = max( 0, tmp1 );
                    tmp2 = min( tmp1, mbloc ); npq += ( tmp2 = min( tmp2, kb ) );
                    mult( *m, *n, ALPHA, aptrd, LDA, B, LDB );
                }
                if( ( kb  -= tmp2 ) == 0 ) return( npq );

                /*
                *  Keep going south until there are no more blocks owning diagonals
                */
                lcmt -= pmb; mblkd--; aptrd += mbloc * inca; B += tmp2 * incb;
            }

            /*
            *  I am done with this column of the LCM table. Go to the next column ...
            */
            lcmt00 += qnb; nblks--;

        /*
        *  ... until there are no more columns.
        */
        } while( nblks > 0 );

        /*
        *  Return the number of diagonals found.
        */
        return( npq );
    }
    else
    {
        /*
        *  (un)packing using columns of virtual matrix
        */
        if( rows )
        {
            /*
            *  (un)packing rows of mn by k array A
            */
            inca = size;
            incb = ( notran ? size : LDB * size );
            m    = &tmp2;
            n    = &K;
        }
        else
        {
            /*
            *  (un)packing columns of k by mn array A
            */
            inca = LDA * size;
            incb = ( notran ? LDB * size : size );
            m    = &K;
            n    = &tmp2;
        }
        kb  = MN;

        /*
        *  From the (un)packing point of view the only valuable shortcut is when the
        *  virtual grid and the blocks are square, and the offset is zero or the grid
        *  is 1x1.
        */
        if( ( ( lcmt00 == 0 ) && ( VM->imb1 == VM->inb1 ) && ( mb == nb ) &&
                ( nprow == npcol ) ) || ( ( nprow == 1 ) && ( npcol == 1 ) ) )
        {
            if(  VM->prow == VM->pcol )
            {
                npq = ( ( nblks <  2 ) ? inbloc :
                        inbloc + ( nblks - 2 ) * nb + lnbloc );
                npq = min( npq, kb );
                if( rows ) mult( npq, K, ALPHA, A, LDA, B, LDB );
                else       mult( K, npq, ALPHA, A, LDA, B, LDB );
            }
            return( npq );
        }
        pmb = nprow * mb;
        qnb = npcol * nb;

        /*
        *  Handle separately the first row and/or column of the LCM table. Update the
        *  LCM value of the curent block lcmt00, as well as the number of rows and
        *  columns mblks and nblks remaining in the LCM table.
        */
        GoSouth = ( lcmt00 > iupp );
        GoEast  = ( lcmt00 < ilow );

        if( !( GoSouth ) && !( GoEast ) )
        {
            /*
            *  The upper left block owns diagonal entries lcmt00 >= ilow && lcmt00 <= iupp
            */
            if( lcmt00 >= 0 )
            {
                tmp1 = imbloc - lcmt00; tmp1 = max( 0, tmp1 );
                tmp2 = min( tmp1, inbloc ); npq += ( tmp2 = min( tmp2, kb ) );
                mult( *m, *n, ALPHA, A, LDA, B, LDB );
            }
            else
            {
                tmp1 = inbloc + lcmt00; tmp1 = max( 0, tmp1 );
                tmp2 = min( tmp1, imbloc ); npq += ( tmp2 = min( tmp2, kb ) );
                mult( *m, *n, ALPHA, A-lcmt00*inca, LDA, B, LDB );
            }
            if( ( kb -= tmp2 ) == 0 ) return( npq );
            B += tmp2 * incb;

            /*
            *  Decide whether one should go south or east in the table: Go east if
            *  the block below the current one only owns lower entries. If this block,
            *  however, owns diagonals, then go south.
            */
            GoSouth = !( GoEast = ( ( lcmt00 - ( iupp - upp + pmb ) ) < ilow ) );
        }

        if( GoSouth )
        {
            /*
            *  Go one step south in the LCM table. Adjust the current LCM value as well as
            *  the pointer to B. The pointer to A remains unchanged.
            */
            lcmt00 -= iupp - upp + pmb; mblks--;

            /*
            *  While there are blocks remaining that own upper entries, keep going south.
            *  Adjust the current LCM value as well as the pointer to B accordingly.
            */
            while( mblks && ( lcmt00 > upp ) ) { lcmt00 -= pmb; mblks--; }

            /*
            *  Return if no more row in the LCM table.
            */
            if( mblks <= 0 ) return( npq );

            /*
            *  lcmt00 <= upp. The current block owns either diagonals or lower entries.
            *  Save the current position in the LCM table. After this column has been
            *  completely taken care of, re-start from this row and the next column of
            *  the LCM table.
            */
            lcmt  = lcmt00; mblkd = mblks;

            while( mblkd && ( lcmt >= ilow ) )
            {
                /*
                *  A block owning diagonals lcmt00 >= ilow && lcmt00 <= upp has been found.
                */
                mbloc = ( ( mblkd == 1 ) ? lmbloc : mb );
                if( lcmt >= 0 )
                {
                    tmp1 = mbloc - lcmt; tmp1 = max( 0, tmp1 );
                    tmp2 = min( tmp1, inbloc ); npq += ( tmp2 = min( tmp2, kb ) );
                    mult( *m, *n, ALPHA, A, LDA, B, LDB );
                }
                else
                {
                    tmp1 = inbloc + lcmt; tmp1 = max( 0, tmp1 );
                    tmp2 = min( tmp1, mbloc ); npq += ( tmp2 = min( tmp2, kb ) );
                    mult( *m, *n, ALPHA, A-lcmt*inca, LDA, B, LDB );
                }
                if( ( kb  -= tmp2 ) == 0 ) return( npq );

                /*
                *  Keep going south until there are no more blocks owning diagonals
                */
                lcmt  -= pmb; mblkd--; B += tmp2 * incb;
            }

            /*
            *  I am done with the first column of the LCM table. Go to the next column.
            */
            lcmt00 += low - ilow + qnb; nblks--; A += inbloc * inca;
        }
        else if( GoEast )
        {
            /*
            *  Go one step east in the LCM table. Adjust the current LCM value as
            *  well as the pointer to A. The pointer to B remains unchanged.
            */
            lcmt00 += low - ilow + qnb; nblks--; A += inbloc * inca;

            /*
            *  While there are blocks remaining that own lower entries, keep going east
            *  in the LCM table. Adjust the current LCM value as well as the pointer to
            *  A accordingly.
            */
            while( nblks && ( lcmt00 < low ) )
            { lcmt00 += qnb; nblks--; A += nb * inca; }

            /*
            *  Return if no more column in the LCM table.
            */
            if( nblks <= 0 ) return( npq );

            /*
            *  lcmt00 >= low. The current block owns either diagonals or upper entries. Save
            *  the current position in the LCM table. After this row has been completely
            *  taken care of, re-start from this column and the next row of the LCM table.
            */
            lcmt  = lcmt00; nblkd = nblks; aptrd = A;

            while( nblkd && ( lcmt <= iupp ) )
            {
                /*
                *  A block owning diagonals lcmt00 >= low && lcmt00 <= iupp has been found.
                */
                nbloc = ( ( nblkd == 1 ) ? lnbloc : nb );
                if( lcmt >= 0 )
                {
                    tmp1 = imbloc - lcmt; tmp1 = max( 0, tmp1 );
                    tmp2 = min( tmp1, nbloc ); npq += ( tmp2 = min( tmp2, kb ) );
                    mult( *m, *n, ALPHA, aptrd, LDA, B, LDB );
                }
                else
                {
                    tmp1 = nbloc + lcmt; tmp1 = max( 0, tmp1 );
                    tmp2 = min( tmp1, imbloc ); npq += ( tmp2 = min( tmp2, kb ) );
                    mult( *m, *n, ALPHA, aptrd-lcmt*inca, LDA, B, LDB );
                }
                if( ( kb  -= tmp2 ) == 0 ) return( npq );

                /*
                *  Keep going east until there are no more blocks owning diagonals.
                */
                lcmt += qnb; nblkd--; aptrd += nbloc * inca; B += tmp2 * incb;
            }

            /*
            *  I am done with the first row of the LCM table. Go to the next row.
            */
            lcmt00 -= iupp - upp + pmb; mblks--;
        }

        /*
        *  Loop over the remaining columns of the LCM table.
        */
        do
        {
            /*
            *  If the current block does not have diagonal elements, find the closest one in
            *  the LCM table having some.
            */
            if( ( lcmt00 < low ) || ( lcmt00 > upp ) )
            {
                while( mblks && nblks )
                {
                    while( mblks && ( lcmt00 > upp ) )
                    { lcmt00 -= pmb; mblks--; }
                    if( lcmt00 >= low ) break;
                    while( nblks && ( lcmt00 < low ) )
                    { lcmt00 += qnb; nblks--; A += nb*inca; }
                    if( lcmt00 <= upp ) break;
                }
            }
            if( !( mblks ) || !( nblks ) ) return( npq );

            /*
            *  The current block owns diagonals. Save the current position in the LCM table.
            *  After this column has been completely taken care of, re-start from this row
            *  and the next column in the LCM table.
            */
            nbloc = ( ( nblks == 1 ) ? lnbloc : nb );
            lcmt = lcmt00; mblkd = mblks;

            while( mblkd && ( lcmt >= low ) )
            {
                /*
                *  A block owning diagonals lcmt00 >= low && lcmt00 <= upp has been found.
                */
                mbloc = ( ( mblkd == 1 ) ? lmbloc : mb );
                if( lcmt >= 0 )
                {
                    tmp1 = mbloc - lcmt; tmp1 = max( 0, tmp1 );
                    tmp2 = min( tmp1, nbloc ); npq += ( tmp2 = min( tmp2, kb ) );
                    mult( *m, *n, ALPHA, A, LDA, B, LDB );
                }
                else
                {
                    tmp1 = nbloc + lcmt; tmp1 = max( 0, tmp1 );
                    tmp2 = min( tmp1, mbloc ); npq += ( tmp2 = min( tmp2, kb ) );
                    mult( *m, *n, ALPHA, A-lcmt*inca, LDA, B, LDB );
                }
                if( ( kb  -= tmp2 ) == 0 ) return( npq );

                /*
                *  Keep going south until there are no more blocks owning diagonals
                */
                lcmt -= pmb; mblkd--; B += tmp2 * incb;
            }

            /*
            *  I am done with this column of the LCM table. Go to the next column ...
            */
            lcmt00 += qnb; nblks--; A += nbloc * inca;

        /*
        *  ... until there are no more columns.
        */
        } while( nblks > 0 );

        /*
        *  Return the number of diagonals found.
        */
        return( npq );
    }

}