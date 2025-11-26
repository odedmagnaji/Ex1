import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Test;

/**
 * JUnit tests for Ex1 polynomial methods (JUnit 5 / Jupiter).
 * Each test checks a specific method of the Ex1 class.
 */
public class Ex1Test {

    /**
     * Tests f(x) on a simple linear polynomial p(x) = 2 + 3x.
     * Verifies that evaluation at a few x values is correct.
     */
    @Test
    public void testF_simple() {
        double[] p = {3, 4}; // p(x) = 3 + 4x
        assertEquals(3.0, Ex1.f(p, 0), Ex1.EPS);
        assertEquals(7.0, Ex1.f(p, 1), Ex1.EPS);  // 3 + 4*1
        assertEquals(11.0, Ex1.f(p, 2), Ex1.EPS);  // 3 + 4*2
    }

    /**
     * Tests the root_rec function on the simplest possible polynomial: f(x) = x.
     * The polynomial is represented as {0, 1} → f(x) = 0 + 1*x.
     *
     * Root of f(x)=x is exactly x=0.
     * Searching in the interval [-1, 1] guarantees a sign change:
     *   f(-1) = -1 < 0
     *   f(1)  =  1 > 0
     *
     * Since the bisection method (root_rec) always checks the midpoint,
     * it will very quickly reach x=0 — where |f(0)| = 0 < EPS.
     *
     * The test verifies that root_rec returns a value very close to 0.
     */
    @Test
    public void testRootRec_linearSimple() {
        double[] p = {0, 1}; // f(x) = x
        double root = Ex1.root_rec(p, -1, 1, Ex1.EPS);
        assertEquals(0.0, root, 0.001);
    }


    /**
     * Tests PolynomFromPoints for the linear case (2 points).
     * Uses points on the line y = 2x + 1 and checks that
     * the returned polynomial is {1,2}.
     */
    @Test
    public void testPolynomFromPoints_linear() {
        double[] xx = {0, 1};
        double[] yy = {1, 3}; // points on y = 2x + 1
        double[] p = Ex1.PolynomFromPoints(xx, yy);
        double[] expected = {1.0, 2.0}; // 1 + 2x
        assertArrayEquals(expected, p, Ex1.EPS);
    }

    /**
     * Tests PolynomFromPoints for the quadratic case (3 points).
     * Uses points on y = x^2 and checks that the returned polynomial is {0,0,1}.
     */
    @Test
    public void testPolynomFromPoints_quadratic() {
        double[] xx = {0, 1, 2};
        double[] yy = {0, 1, 4}; // points on y = x^2
        double[] p = Ex1.PolynomFromPoints(xx, yy);
        double[] expected = {0.0, 0.0, 1.0}; // x^2
        assertArrayEquals(expected, p, 1e-3);
    }

    /**
     * Tests that two identical polynomials are considered equal.
     */
    @Test
    public void testEquals_basic() {
        double[] p1 = {2.0, 3.0, 0.0};
        double[] p2 = {2.0, 3.0, 0.0};
        assertTrue(Ex1.equals(p1, p2));
    }
    /**
     * Tests that trailing zeros do not affect polynomial equality.
     */
    @Test
    public void testEquals_trailingZeros() {
        double[] p1 = {2.0, 3.0, 0.0};
        double[] p2 = {2.0, 3.0};
        assertTrue(Ex1.equals(p1, p2));
    }
    /**
     * Tests that polynomials with different coefficients are not considered equal.
     */
    @Test
    public void testEquals_differentCoefficients() {
        double[] p1 = {2.0, 3.0};
        double[] p2 = {2.0, 3.1};
        assertFalse(Ex1.equals(p1, p2));
    }



    /**
     * Tests that poly(...) and getPolynomFromString(...) are inverse operations.
     * Converts an array to String and back and checks equality.
     */
    @Test
    public void testPoly_and_getPolynomFromString_inverse() {
        double[] p = {2.0, -3.0, 1.0}; // 1x^2 -3x +2
        String s = Ex1.poly(p);
        double[] back = Ex1.getPolynomFromString(s);
        assertTrue(Ex1.equals(p, back));
    }
    /**
     * Tests that poly(...) produces the expected String representation
     * for a general polynomial with positive and negative coefficients,
     * skipping zero terms.
     */
    @Test
    public void testPoly_format() {
        double[] p = {2.0, 0.0, 3.1, -1.2}; // -1.2x^3 +3.1x^2 +2.0
        String s = Ex1.poly(p);
        assertEquals("-1.2x^3 +3.1x^2 +2.0", s);
    }
    /**
     * Tests that the zero polynomial is represented as "0".
     */
    @Test
    public void testPoly_zero() {
        double[] p = {0.0, 0.0, 0.0};
        String s = Ex1.poly(p);
        assertEquals("0", s);
    }


    /**
     * Tests sameValue(p1,p2,...) on p1(x)=x and p2(x)=0 in [-1,1].
     * The intersection point should be close to x=0.
     */
    @Test
    public void testSameValue_simple() {
        double[] p1 = {0.0, 1.0}; // y = x
        double[] p2 = {0.0};      // y = 0
        double x = Ex1.sameValue(p1, p2, -1.0, 1.0, Ex1.EPS);
        assertEquals(0.0, x, 0.01);
    }

    /**
     * Tests length(p,x1,x2,segments) on the zero function p(x)=0.
     * The length from x=0 to x=10 should be exactly 10.
     */
    @Test
    public void testLength_zeroFunction() {
        double[] p = {0.0}; // y = 0
        double len = Ex1.length(p, 0.0, 10.0, 100);
        assertEquals(10.0, len, 1e-3);
    }
    /**
     * Tests length with only one segment: we get direct straight-line distance.
     */
    @Test
    public void testLength_oneSegment() {
        double[] p = {0.0, 1.0}; // y = x
        double len = Ex1.length(p, 0.0, 3.0, 1);
        // One straight segment between (0,0) and (3,3)
        assertEquals(Math.sqrt(18), len, 1e-6);
    }

    /**
     * Tests area(p1,p2,...) for p1(x)=x, p2(x)=0 on [0,1].
     * The exact area is 1/2, approximated using many trapezoids.
     */
    @Test
    public void testArea_simpleTriangle() {
        double[] p1 = {0.0, 1.0}; // y = x
        double[] p2 = {0.0};      // y = 0
        double a = Ex1.area(p1, p2, 0.0, 1.0, 1000);
        assertEquals(0.5, a, 0.01);
    }
    /**
     * Tests that the area between two identical polynomials is zero (up to a small tolerance).
     */
    @Test
    public void testArea_identicalPolynomials() {
        double[] p1 = {1.0, -2.0, 3.0}; // some quadratic
        double[] p2 = {1.0, -2.0, 3.0}; // exactly the same
        double a = Ex1.area(p1, p2, -5.0, 5.0, 1000);
        assertEquals(0.0, a, 1e-3);
    }

    /**
     * Tests getPolynomFromString(...) on the example from the assignment:
     * "-1.0x^2 +3.0x +2.0" -> {2.0, 3.0, -1.0}.
     */
    @Test
    public void testGetPolynomFromString_example() {
        String s = "-1.0x^2 +3.0x +2.0";
        double[] p = Ex1.getPolynomFromString(s);
        double[] expected = {2.0, 3.0, -1.0};
        assertArrayEquals(expected, p, Ex1.EPS);
    }
    /**
     * Tests parsing of a polynomial without spaces: "3x^2-5x+1".
     */
    @Test
    public void testGetPolynomFromString_noSpaces() {
        String s = "3x^2-5x+1";
        double[] p = Ex1.getPolynomFromString(s);
        double[] expected = {1.0, -5.0, 3.0};
        assertArrayEquals(expected, p, Ex1.EPS);
    }
    /**
     * Tests parsing when coefficients are implicit (e.g., "x^2-x+1").
     */
    @Test
    public void testGetPolynomFromString_implicitCoefficients() {
        String s = "x^2 - x + 1";
        double[] p = Ex1.getPolynomFromString(s);
        double[] expected = {1.0, -1.0, 1.0};
        assertArrayEquals(expected, p, Ex1.EPS);
    }
    /**
     * Tests parsing of a constant polynomial: "5" -> {5.0}.
     */
    @Test
    public void testGetPolynomFromString_constant() {
        String s = "5";
        double[] p = Ex1.getPolynomFromString(s);
        double[] expected = {5.0};
        assertArrayEquals(expected, p, Ex1.EPS);
    }
    /**
     * Tests parsing of zero polynomial: "0" -> {0}.
     */
    @Test
    public void testGetPolynomFromString_zero() {
        assertArrayEquals(Ex1.ZERO, Ex1.getPolynomFromString("0"), Ex1.EPS);
        assertArrayEquals(Ex1.ZERO, Ex1.getPolynomFromString(""), Ex1.EPS);
        assertArrayEquals(Ex1.ZERO, Ex1.getPolynomFromString("   "), Ex1.EPS);
    }
    /**
     * Tests parsing of a polynomial with a high degree (e.g., x^7 - 2x^3 + 1).
     */
    @Test
    public void testGetPolynomFromString_highPower() {
        String s = "x^7 - 2x^3 + 1";
        double[] p = Ex1.getPolynomFromString(s);
        double[] expected = {1.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, 1.0};
        assertArrayEquals(expected, p, Ex1.EPS);
    }
    /**
     * Tests parsing of alternating signs: "-2x + 4x^3 - 5".
     */
    @Test
    public void testGetPolynomFromString_alternatingSigns() {
        String s = "-2x + 4x^3 - 5";
        double[] p = Ex1.getPolynomFromString(s);
        double[] expected = {-5.0, -2.0, 0.0, 4.0};
        assertArrayEquals(expected, p, Ex1.EPS);
    }

    /**
     * Tests add(p1,p2) on two simple polynomials.
     * (1 + 2x) + (3 - x) = 4 + 1x.
     */
    @Test
    public void testAdd_simple() {
        double[] p1 = {1.0, 2.0};   // 1 + 2x
        double[] p2 = {3.0, -1.0};  // 3 - x
        double[] sum = Ex1.add(p1, p2);
        double[] expected = {4.0, 1.0}; // 4 + 1x
        assertArrayEquals(expected, sum, Ex1.EPS);
    }
    /**
     * Tests add(p1,p2) when polynomials have different lengths.
     * (2 + 3x) + (1 + x + x^2) = 3 + 4x + x^2.
     */
    @Test
    public void testAdd_differentLengths() {
        double[] p1 = {2.0, 3.0};      // 2 + 3x
        double[] p2 = {1.0, 1.0, 1.0}; // 1 + x + x^2
        double[] sum = Ex1.add(p1, p2);
        double[] expected = {3.0, 4.0, 1.0};
        assertArrayEquals(expected, sum, Ex1.EPS);
    }
    /**
     * Tests that adding the zero polynomial does not change the polynomial.
     */
    @Test
    public void testAdd_withZero() {
        double[] p1 = {5.0, -1.0, 2.0};
        double[] zero = Ex1.ZERO;
        assertArrayEquals(p1, Ex1.add(p1, zero), Ex1.EPS);
        assertArrayEquals(p1, Ex1.add(zero, p1), Ex1.EPS);
    }

    /**
     * Tests mul(p1,p2) on (1 + x)*(1 + x).
     * Expected result is (1 + 2x + x^2).
     */
    @Test
    public void testMul_simple() {
        double[] p1 = {1.0, 1.0};  // 1 + x
        double[] p2 = {1.0, 1.0};  // 1 + x
        double[] prod = Ex1.mul(p1, p2);
        double[] expected = {1.0, 2.0, 1.0}; // 1 + 2x + x^2
        assertArrayEquals(expected, prod, Ex1.EPS);
    }
    /**
     * Multiplying any polynomial by ZERO must return ZERO.
     */
    @Test
    public void testMul_zero() {
        double[] p1 = {1.0, 2.0, 3.0};
        double[] zero = {0.0};

        assertArrayEquals(zero, Ex1.mul(p1, zero), Ex1.EPS);
        assertArrayEquals(zero, Ex1.mul(zero, p1), Ex1.EPS);
    }
    /**
     * Multiplying a polynomial by a constant scalar.
     * (2) * (1 + x + x^2) = 2 + 2x + 2x^2
     */
    @Test
    public void testMul_constant() {
        double[] constant = {2.0};
        double[] poly = {1.0, 1.0, 1.0};
        double[] expected = {2.0, 2.0, 2.0};

        assertArrayEquals(expected, Ex1.mul(constant, poly), Ex1.EPS);
    }
    /**
     * Multiplying polynomials of different lengths.
     * (1 + x + x^2) * (3 + 0x) = 3 + 3x + 3x^2
     */
    @Test
    public void testMul_differentLengths() {
        double[] p1 = {1.0, 1.0, 1.0};
        double[] p2 = {3.0};
        double[] expected = {3.0, 3.0, 3.0};

        assertArrayEquals(expected, Ex1.mul(p1, p2), Ex1.EPS);
    }
    /**
     * Multiplying by 1 should return the same polynomial.
     */
    @Test
    public void testMul_identity() {
        double[] one = {1.0};
        double[] p = {5.0, -3.0, 2.0};

        assertArrayEquals(p, Ex1.mul(one, p), Ex1.EPS);
    }

    /**
     * Tests derivative(...) on a non-constant polynomial:
     * p(x)=2+3x+4x^2 -> p'(x)=3+8x.
     */
    @Test
    public void testDerivative_simple() {
        double[] p = {2.0, 3.0, 4.0};
        double[] dp = Ex1.derivative(p);
        double[] expected = {3.0, 8.0};
        assertArrayEquals(expected, dp, Ex1.EPS);
    }

    /**
     * Tests derivative(...) on a constant polynomial.
     * The derivative of p(x)=5 should be the zero polynomial.
     */
    @Test
    public void testDerivative_constant() {
        double[] p = {5.0}; // constant polynomial
        double[] dp = Ex1.derivative(p);
        assertArrayEquals(Ex1.ZERO, dp, Ex1.EPS);
    }
}
