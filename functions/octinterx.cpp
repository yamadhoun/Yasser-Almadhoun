#include <octave/oct.h>
#include <math.h>


static bool any_bad_argument(const octave_value_list& args)
{
    if (!args(0).is_real_matrix())
    {
        error("octinterx: expecting x1 (arg 1) to be a real matrix");
        return true;
    }

    if (!args(1).is_real_matrix())
    {
        error("octinterx: expecting y1 (arg 1) to be a real matrix");
        return true;
    }

    if (args.length()==4){
       if (!args(2).is_real_matrix())
       {
           error("octinterx: expecting x2 (arg 1) to be a real matrix");
           return true;
       }

       if (!args(3).is_real_matrix())
       {
           error("octinterx: expecting y2 (arg 1) to be a real matrix");
           return true;
       }
    }

    return false;
}




DEFUN_DLD (octinterx, args, nargout,
           "[x0,y0,iout,jout] = octinterx(x1,y1,x2,y2)\n\
           \n\
           Octave C++ implementation for intersections.m\n\
           Computes the (x,y) locations where two curves intersect.  The curves\n\
           can be broken with NaNs or have vertical segments.\n\
           \n\
           where X1 and Y1 are equal-length vectors of at least two points and\n\
           represent curve 1.  Similarly, X2 and Y2 represent curve 2.\n\
           X0 and Y0 are column vectors containing the points at which the two\n\
           curves intersect.\n\
           For each element of the vector I, I(k) = (segment number of (X1,Y1)) +\n\
           (how far along this segment the intersection is).  For example, if I(k) =\n\
           45.25 then the intersection lies a quarter of the way between the line\n\
           segment connecting (X1(45),Y1(45)) and (X1(46),Y1(46)).  Similarly for\n\
           the vector J and the segments in (X2,Y2).\n\
           \n\
           You can also get intersections of a curve with itself.  Simply pass in\n\
           only one curve")
{
        // input check
        if (args.length() != 4 && args.length() != 2)
               print_usage();

        if ( any_bad_argument(args))
               print_usage();

        // Read input coordinates
        NDArray x1 = args(0).array_value();
        NDArray y1 = args(1).array_value();
        NDArray x2 = args(2).array_value();
        NDArray y2 = args(3).array_value();

        // Find intersections
        float eps=1e-4;

        // check circularity
        octave_idx_type np = x1.columns()-1;
        if (abs(x1(0)-x1(np))<eps && abs(y1(0)-y1(np))<eps) {
            int a = 0;


           }


        // Return empty matrices for any outputs
        octave_value_list retval (nargout);

        return retval;
}


