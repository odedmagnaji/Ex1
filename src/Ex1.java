/**
 * Introduction to Computer Science 2026, Ariel University,
 * Ex1: arrays, static functions and JUnit
 *
 * This class represents a set of static methods on polynomial functions - represented as an array of doubles.
 * The array {0.1, 0, -3, 0.2} represents the following polynomial function: 0.2x^3-3x^2+0.1
 * This is the main Class you should implement.
 *
 * @author boaz.benmoshe
 */
public class Ex1 {

    /** Epsilon value for numerical computation, it serves as a "close enough" threshold. */
    public static final double EPS = 0.001; // the epsilon to be used for the root approximation.
    /** The zero polynomial function is represented as an array with a single (0) entry. */
    public static final double[] ZERO = {0};

    /**
     * Computes the f(x) value of the polynomial function at x.
     * @param poly - polynomial function
     * @param x - the x value
     * @return f(x) - the polynomial function value at x.
     */
    public static double f(double[] poly, double x) {
        // Evaluate polynomial: sum_i poly[i] * x^i
        double ans = 0;
        for (int i = 0; i < poly.length; i++) {
            double c = Math.pow(x, i);
            ans += c * poly[i];
        }
        return ans;
    }

    /**
     * Given a polynomial function (p), a range [x1,x2] and an epsilon eps.
     * This function computes an x value (x1<=x<=x2) for which |p(x)| < eps,
     * assuming p(x1)*p(x2) <= 0.
     * Implemented recursively using bisection method.
     */
    public static double root_rec(double[] p, double x1, double x2, double eps) {
        // Recursive bisection to find root of p(x) in [x1,x2]
        double f1 = f(p, x1);
        double x12 = (x1 + x2) / 2;
        double f12 = f(p, x12);
        if (Math.abs(f12) < eps) {
            return x12;
        }
        if (f12 * f1 <= 0) {
            return root_rec(p, x1, x12, eps);
        } else {
            return root_rec(p, x12, x2, eps);
        }
    }
    /**
     * Builds a polynomial from 2 or 3 points on its graph.
     * Representation: the returned array is {c, b, a} so that:
     *      p(x) = c + b*x + a*x^2
     *
     * Supported cases:
     *  - 2 points: returns a linear polynomial (degree 1).
     *  - 3 points: returns a quadratic polynomial (degree 2),
     *              using closed-form formulas (based on Lagrange / StackOverflow).
     *
     * If the input is invalid (null, different lengths, or length not 2 or 3),
     * the function returns null.
     */
    public static double[] PolynomFromPoints(double[] xx, double[] yy) {
        double[] ans = null;
        int lx = xx.length;
        int ly = yy.length;

        // basic validation: same length, not null, and 2 or 3 points only
        if (xx != null && yy != null && lx == ly && lx > 1 && lx < 4) {

            // ----- Linear case: 2 points → y = m*x + b -----
            if (lx == 2) {
                double x1 = xx[0], x2 = xx[1];
                double y1 = yy[0], y2 = yy[1];

                // slope m = (y2 - y1) / (x2 - x1)
                double m = (y2 - y1) / (x2 - x1);
                // intercept b = y1 - m*x1
                double b = y1 - m * x1;

                // polynomial representation: {b, m} → p(x) = b + m*x
                ans = new double[2];
                ans[0] = b;
                ans[1] = m;
            }

            // ----- Quadratic case: 3 points → y = a*x^2 + b*x + c -----
            else if (lx == 3) {
                double x1 = xx[0], x2 = xx[1], x3 = xx[2];
                double y1 = yy[0], y2 = yy[1], y3 = yy[2];

                // denominator (non-zero if x1, x2, x3 are distinct)
                double denom = (x1 - x2) * (x1 - x3) * (x2 - x3);

                // If denom is ~0, points are degenerate → no unique quadratic
                if (Math.abs(denom) < EPS) {
                    return null;
                }

                // Closed-form formulas for a, b, c (from 3-point quadratic interpolation)
                double a = ( x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2) ) / denom;
                double b = ( x3 * x3 * (y1 - y2)
                        + x2 * x2 * (y3 - y1)
                        + x1 * x1 * (y2 - y3) ) / denom;
                double c = ( x2 * x3 * (x2 - x3) * y1
                        + x3 * x1 * (x3 - x1) * y2
                        + x1 * x2 * (x1 - x2) * y3 ) / denom;

                // polynomial representation: {c, b, a} → p(x) = c + b*x + a*x^2
                ans = new double[3];
                ans[0] = c;
                ans[1] = b;
                ans[2] = a;
            }
        }
        return ans;
    }


    /**
     * Two polynomials functions are equal if and only if they have the same coefficients
     * up to EPS, ignoring trailing zeros.
     * @param p1 first polynomial function
     * @param p2 second polynomial function
     * @return true iff p1 represents the same polynomial function as p2.
     */
    public static boolean equals(double[] p1, double[] p2) {
        // Compare coefficient by coefficient, allowing small EPS difference
        int maxLen = Math.max(p1.length, p2.length);
        for (int i = 0; i < maxLen; i++) {
            double c1 = (i < p1.length) ? p1[i] : 0;
            double c2 = (i < p2.length) ? p2[i] : 0;
            if (Math.abs(c1 - c2) > EPS) {
                return false;
            }
        }
        return true;
    }

    /**
     * Computes a String representing the polynomial function.
     * For example the array {2,0,3.1,-1.2} will be presented as "-1.2x^3 +3.1x^2 +2.0"
     * @param poly the polynomial function represented as an array of doubles
     * @return String representing the polynomial function.
     */
    public static String poly(double[] poly) {
        // Build a human-readable string from highest to lowest degree
        String ans = "";
        if (poly == null || poly.length == 0) {
            return "0";
        }
        boolean first = true;
        for (int i = poly.length - 1; i >= 0; i--) {
            double c = poly[i];
            if (Math.abs(c) < EPS) {
                continue; // skip near-zero coefficients
            }
            String term = "";
            String sign = (c < 0) ? "-" : "+";
            double absC = Math.abs(c);

            if (first) {
                // First non-zero term: only minus if needed
                if (c < 0) {
                    term += "-";
                }
            } else {
                // Next terms: include sign with space
                term += " " + sign;
            }

            if (i == 0) {
                // Constant term
                term += absC;
            } else if (i == 1) {
                // Linear term
                term += absC + "x";
            } else {
                // Higher degree term
                term += absC + "x^" + i;
            }
            ans += term;
            first = false;
        }
        if (first) {
            // All coefficients were zero
            return "0";
        }
        return ans;
    }

    /**
     * Given two polynomial functions (p1,p2), a range [x1,x2] and an epsilon eps.
     * This function computes an x value (x1<=x<=x2) for which |p1(x) -p2(x)| < eps,
     * assuming (p1(x1)-p2(x1)) * (p1(x2)-p2(x2)) <= 0.
     * Uses recursive bisection on the function g(x) = p1(x)-p2(x).
     */
    public static double sameValue(double[] p1, double[] p2, double x1, double x2, double eps) {
        // Bisection on g(x) = f(p1,x) - f(p2,x)
        double g1 = f(p1, x1) - f(p2, x1);
        double mid = (x1 + x2) / 2.0;
        double gm = f(p1, mid) - f(p2, mid);

        if (Math.abs(gm) < eps) {
            return mid;
        }
        if (g1 * gm <= 0) {
            return sameValue(p1, p2, x1, mid, eps);
        } else {
            return sameValue(p1, p2, mid, x2, eps);
        }
    }

    /**
     * Given a polynomial function (p), a range [x1,x2] and an integer with the number (n) of sample points.
     * This function computes an approximation of the length of the function between f(x1) and f(x2)
     * using n segments and computing the path length between neighboring points.
     * Implemented iteratively (non-recursive).
     */
    public static double length(double[] p, double x1, double x2, int numberOfSegments) {
        // Approximate arc length using piecewise linear segments
        double h = (x2 - x1) / numberOfSegments;
        double prevX = x1;
        double prevY = f(p, x1);
        double sum = 0.0;

        for (int i = 1; i <= numberOfSegments; i++) {
            double x = x1 + i * h;
            double y = f(p, x);
            double dx = x - prevX;
            double dy = y - prevY;
            sum += Math.sqrt(dx * dx + dy * dy);
            prevX = x;
            prevY = y;
        }
        return sum;
    }

    /**
     * Given two polynomial functions (p1,p2), a range [x1,x2] and an integer representing
     * the number of Trapezoids between the functions.
     * This function computes an approximation of the area between the polynomial functions
     * within the x-range using the trapezoidal rule on |p1(x)-p2(x)|.
     */
    public static double area(double[] p1, double[] p2, double x1, double x2, int numberOfTrapezoid) {
        // Approximate area between graphs: integral of |p1(x)-p2(x)| using trapezoids
        if (numberOfTrapezoid <= 0 || x1 == x2) {
            return 0;
        }
        double h = (x2 - x1) / numberOfTrapezoid;
        double sum = 0.0;

        for (int i = 0; i < numberOfTrapezoid; i++) {
            double xLeft = x1 + i * h;
            double xRight = xLeft + h;

            double yLeft = Math.abs(f(p1, xLeft) - f(p2, xLeft));
            double yRight = Math.abs(f(p1, xRight) - f(p2, xRight));

            sum += (yLeft + yRight) * 0.5 * h;
        }
        return sum;
    }

    /**
     * This function computes the array representation of a polynomial function from a String
     * representation. Example:
     * "-1.0x^2 +3.0x +2.0"  -->  {2.0, 3.0, -1.0}
     * @param p - a String representing polynomial function.
     * @return double[] representation of the polynomial.
     */
    public static double[] getPolynomFromString(String p) {
        // Parse polynomial string into coefficient array
        double[] ans = ZERO;
        if (p == null) {
            return ZERO;
        }
        String s = p.replace(" ", "");
        if (s.length() == 0 || s.equals("0")) {
            return ZERO;
        }
        // Ensure first char is a sign
        if (s.charAt(0) != '+' && s.charAt(0) != '-') {
            s = "+" + s;
        }
        java.util.ArrayList<String> terms = new java.util.ArrayList<String>();
        StringBuilder curr = new StringBuilder();
        for (int i = 0; i < s.length(); i++) {
            char c = s.charAt(i);
            if ((c == '+' || c == '-') && i != 0) {
                terms.add(curr.toString());
                curr = new StringBuilder();
                curr.append(c);
            } else {
                curr.append(c);
            }
        }
        if (curr.length() > 0) {
            terms.add(curr.toString());
        }

        // First pass: find max power
        int maxPow = 0;
        for (String term : terms) {
            if (term.length() == 0) continue;
            String rest = term.substring(1); // without sign
            int pow;
            if (rest.contains("x")) {
                if (rest.contains("^")) {
                    String powStr = rest.substring(rest.indexOf("^") + 1);
                    pow = Integer.parseInt(powStr);
                } else {
                    pow = 1;
                }
            } else {
                pow = 0;
            }
            if (pow > maxPow) {
                maxPow = pow;
            }
        }

        ans = new double[maxPow + 1];
        // Second pass: fill coefficients
        for (String term : terms) {
            if (term.length() == 0) continue;
            int sign = (term.charAt(0) == '-') ? -1 : 1;
            String rest = term.substring(1);
            double coeff;
            int pow;
            if (rest.contains("x")) {
                int idxX = rest.indexOf('x');
                String coeffStr = rest.substring(0, idxX);
                if (coeffStr.equals("")) {
                    coeffStr = "1";
                }
                coeff = Double.parseDouble(coeffStr) * sign;
                if (rest.contains("^")) {
                    String powStr = rest.substring(rest.indexOf("^") + 1);
                    pow = Integer.parseInt(powStr);
                } else {
                    pow = 1;
                }
            } else {
                coeff = Double.parseDouble(rest) * sign;
                pow = 0;
            }
            ans[pow] += coeff;
        }

        // Trim trailing zeros
        int last = ans.length - 1;
        while (last > 0 && Math.abs(ans[last]) < EPS) {
            last--;
        }
        if (last < ans.length - 1) {
            double[] trimmed = new double[last + 1];
            for (int i = 0; i <= last; i++) {
                trimmed[i] = ans[i];
            }
            ans = trimmed;
        }
        if (ans.length == 1 && Math.abs(ans[0]) < EPS) {
            ans = ZERO;
        }
        return ans;
    }

    /**
     * This function computes the polynomial function which is the sum of two polynomial functions (p1,p2)
     * @param p1 first polynomial
     * @param p2 second polynomial
     * @return p1 + p2
     */
    public static double[] add(double[] p1, double[] p2) {
        // Add polynomials coefficient-wise
        int maxLen = Math.max(p1.length, p2.length);
        double[] ans = new double[maxLen];
        for (int i = 0; i < maxLen; i++) {
            double c1 = (i < p1.length) ? p1[i] : 0;
            double c2 = (i < p2.length) ? p2[i] : 0;
            ans[i] = c1 + c2;
        }
        return ans;
    }

    /**
     * This function computes the polynomial function which is the multiplication of two polynoms (p1,p2)
     * @param p1 first polynomial
     * @param p2 second polynomial
     * @return p1 * p2
     */
    public static double[] mul(double[] p1, double[] p2) {
        if ((p1.length == 1 && Math.abs(p1[0]) < EPS) ||
                (p2.length == 1 && Math.abs(p2[0]) < EPS)) {
            return ZERO;
        }
        double[] ans = new double[p1.length + p2.length - 1];
        for (int i = 0; i < p1.length; i++) {
            for (int j = 0; j < p2.length; j++) {
                ans[i + j] += p1[i] * p2[j];
            }
        }
        return ans;
    }

    /**
     * This function computes the derivative of the p0 polynomial function.
     * @param po original polynomial
     * @return derivative polynomial
     */
    public static double[] derivative(double[] po) {
        // Compute derivative: if p(x)=a0+a1x+a2x^2+..., then p'(x)=a1+2a2x+3a3x^2+...
        if (po.length <= 1) {
            return ZERO;
        }
        double[] ans = new double[po.length - 1];
        for (int i = 1; i < po.length; i++) {
            ans[i - 1] = i * po[i];
        }
        return ans;
    }
}
