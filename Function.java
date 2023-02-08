package thesis.CEC2017;                // Change to be apart of relevant package

// File reading
import java.io.File;                   // Import the File class
import java.io.FileNotFoundException;  // Import this class to handle errors
import java.util.Scanner;              // Import the Scanner class to read text files

/**
 * Base class for all CEC2017 benchmark problems.
 * 
 * @author Liam McDevitt
 * 
 * Function implementations are based on implementations by 
 * Noor Awad in C++.
 * 
 * Works for dimensions 30, 50, & 100.
 * The data file names are adjusted to match those defined in the report.
 * 
 * CEC2017 GitHub : https://github.com/P-N-Suganthan/CEC2017-BoundContrained
 * 
 * Date: 2022/05/24 - 2022/05/30
 */

public class Function {

    private int    n;              // Number of dimensions
    private double min;            // Minimum bound in search space 
    private double max;            // Maximum bound in search space
    private double range;          // Range of search space values

    private int funcNum;          // Which function from the suite to define
    private int compNum;          // Number of composition functions

    private double [][] rotation; // Rotation matrix
    private double [][]   oShift; // Optimum shift
    private int    [][] shuffle;  // Shuffle

    private double [] y, z;       // Temporary decision vectors

    /**
     * Initializes a problem (i.e. an objective function to optimize).
     * @param n        Number of dimensions
     * @param funcNum  Function this object is representing
     */
    public Function (int n, int funcNum) {

        // Initializing global variables
        this.n       = n;
        this.funcNum = funcNum;
        this.compNum = 10;
        this.min     = -100.0;
        this.max     = 100.0;
        this.range   = Math.abs(max - min);

        // Load appropriate shift, rotate, and shuffle data depending on the function
        loadData();

    } // Constructor

    /**
     * Loads the appropriate rotation, shift, and/or shuffle data depending on the function.
     */
    private void loadData () {

        switch (funcNum) {
            case 1: // Shifted and Rotated Bent Cigar Function
            case 2: // Shifted and Rotated Zakharov Function
            case 3: // Shifted and Rotated Rosenbrock’s Function
            case 4: // Shifted and Rotated Rastrigin’s Function
            case 5: // Shifted and Rotated Schaffer F7 Function
            case 6: // Shifted and Rotated Lunacek Bi-Rastrigin’s Function
            case 7: // Non-continuous Rotated Rastrigin’s Function
            case 8: // Levy Function
            case 9: // Modified Schwefel’s Function
            case 20: // Composition function 1
            case 21: // Composition function 2
            case 22: // Composition function 3
            case 23: // Composition function 4
            case 24: // Composition function 5
            case 25: // Composition function 6
            case 26: // Composition function 7
            case 27: // Composition function 8
                loadRotationData();
                loadShiftData();
                break;
            case 10: // Hybrid function 1
            case 11: // Hybrid function 2
            case 12: // Hybrid function 3
            case 13: // Hybrid function 4
            case 14: // Hybrid function 5
            case 15: // Hybrid function 6
            case 16: // Hybrid function 7
            case 17: // Hybrid function 8
            case 18: // Hybrid function 9
            case 19: // Hybrid function 10
            case 28: // Composition function 9
            case 29: // Composition function 10
                loadRotationData();
                loadShiftData();
                loadShuffleData();
                break;
        }

    } // loadData

    /**
     * Load rotation matrix data from file.
     * 
     * TODO: Change to relative file path
     */
    private void loadRotationData () {

        // Try to populate the rotation matrix with rotation data file
        try {
            // Initialize rotation matrix
            if (funcNum <= 19) rotation = new double [n][n];
            else               rotation = new double [compNum * n][n];
            // Initialize file
            File file = new File ("src/main/java/thesis/CEC2017/input_data/M_" + funcNum + "_D" + n + ".txt");
            // Initialize file reader
            Scanner fileReader = new Scanner(file);
            // Fill rotation matrix with file data
            for (int i = 0; i < rotation.length; i++)
                for (int j = 0; j < rotation[i].length; j++)
                    rotation[i][j] = fileReader.nextDouble();
            // Close file
            fileReader.close();
        }
        // Catch file not found exception
        catch (FileNotFoundException fnfe) {
            System.out.println("File not found exception.");
        }

    } // loadRotationData

    /**
     * Load shift matrix data from file.
     * 
     * TODO: Change to relative file path
     */
    private void loadShiftData () {

        // Try to populate the rotation matrix with rotation data file
        try {
            // Initialize shift matrix
            if (funcNum <= 19) oShift = new double [1][n];
            else               oShift = new double [compNum][n];
            // Initialize file
            File file = new File ("src/main/java/thesis/CEC2017/input_data/shift_data_" + funcNum + ".txt");
            // Initialize file reader
            Scanner fileReader = new Scanner(file);
            // Fill rotation matrix with file data
            for (int i = 0; i < oShift.length; i++) {
                for (int j = 0; j < oShift[i].length; j++) {
                    oShift[i][j] = fileReader.nextDouble();
                }
                fileReader.nextLine();
            }
            // Close file
            fileReader.close();
        }
        // Catch file not found exception
        catch (FileNotFoundException fnfe) {
            System.out.println("File not found exception.");
        }

    } // loadShiftData

    /**
     * Load shuffle matrix data from file.
     * 
     * TODO: Change to relative file path
     */
    private void loadShuffleData () {

        // Try to populate the rotation matrix with rotation data file
        try {
            // Initialize shuffle matrix
            if      (funcNum >= 10 && funcNum <= 19) shuffle = new int [1][n];
            else if (funcNum == 28 || funcNum == 29) shuffle = new int [compNum][n];
            // Initialize file
            File file = new File ("src/main/java/thesis/CEC2017/input_data/shuffle_data_" + funcNum + "_D" + n + ".txt");
            // Initialize file reader
            Scanner fileReader = new Scanner(file);
            // Fill rotation matrix with file data
            for (int i = 0; i < shuffle.length; i++) {
                for (int j = 0; j < shuffle[i].length; j++) {
                    shuffle[i][j] = fileReader.nextInt();
                }
            }
            // Close file
            fileReader.close();
        }
        // Catch file not found exception
        catch (FileNotFoundException fnfe) {
            System.out.println("File not found exception.");
        }

    } // loadShuffleData

    /**
     * Shifts function's global optimum.
     * @param x       Vector to shift
     * @param xShift  Shifted vector
     * @param funcPos Function's position in transformation matrices
     */
    private void shiftFunc (double [] x, double [] xShift, int funcPos) {
        for (int i = 0; i < x.length; i++)
            xShift[i] = x[i] - oShift[funcPos][i];
    } // shiftFunc

    /**
     * Rotates function.
     * @param x       Vector to rotate
     * @param xShift  Rotated vector
     * @param funcPos Function's position in transformation matrices
     */
    private void rotateFunc (double [] x, double [] xRotate, int funcPos) {
        for (int i = 0; i < x.length; i++) {
            xRotate[i] = 0;
            for (int j = 0; j < x.length; j++) {
                xRotate[i] += x[j] * rotation[i + (x.length * funcPos)][j]; 
            }
        }
    } // rotateFunc

    /**
     * Shift and rotate function.
     * @param x            Vector to shift and rotate
     * @param xShiftRotate Shifted and rotated vector
     * @param shiftRate    Shift rate
     * @param s_flag       Function shift flag
     * @param r_flag       Function rotation flag
     * @param funcPos      Function's position in transformation matrices
     */
    private void shiftRotateFunc (double [] x, double [] xShiftRotate, double shiftRate, int s_flag, int r_flag, int funcPos) {
        if (s_flag == 1) {
            if (r_flag == 1) {
                shiftFunc(x, y, funcPos);
                for (int i = 0; i < x.length; i++) // "Shrink to the original search range"
                y[i] *= shiftRate;
                rotateFunc(y, xShiftRotate, funcPos);
            }
            else {
                shiftFunc(x, xShiftRotate, funcPos);
                for (int i = 0; i < x.length; i++)
                    xShiftRotate[i] *= shiftRate;
            }
        }
        else {
            if (r_flag == 1) {
                for (int i = 0; i < x.length; i++)
                y[i] = x[i] * shiftRate;
                rotateFunc(y, xShiftRotate, funcPos);
            }
            else {
                for (int i = 0; i < x.length; i++)
                    xShiftRotate[i] = x[i] * shiftRate;
            }
        }
    } // shiftRotateFunc

    /**
     * Calculates the fitness of a composition function.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param delta   Delta values
     * @param bias    Bias values
     * @param fit     Fitness values from each function within the composition
     * @param cf_num  Number of functions in the composition.
     * @return        double: objective fitness
     */
    private double cf_cal (double [] x, int nx, double [] delta, double [] bias, double [] fit, int cf_num) {
        int i,j;
        double [] w = new double [cf_num];
        double w_max = 0, w_sum = 0, INF = 1.0e99;
        for (i = 0; i < cf_num; i++) {
            fit[i] += bias[i];
            w[i] = 0;
            for (j = 0; j < nx; j++) {
                w[i] += Math.pow(x[j] - oShift[i][j], 2.0);
            }
            if (w[i] != 0)
                w[i] = Math.pow(1.0 / w[i], 0.5) * Math.exp(-w[i] / 2.0 / nx / Math.pow(delta[i], 2.0));
            else
                w[i] = INF;
            if (w[i] > w_max)
                w_max = w[i];
        }
        for (i = 0; i < cf_num; i++) {
            w_sum = w_sum + w[i];
        }
        if(w_max == 0) {
            for (i = 0; i < cf_num; i++)
                w[i] = 1;
            w_sum = cf_num;
        }
        double result = 0.0;
        for (i = 0; i < cf_num; i++) {
            result = result + w[i] / w_sum * fit[i];
        }
        return result;
    } // cf_cal

    /**
     * High Conditioned Elliptic Function.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param s_flag  Function shift flag
     * @param r_flag  Function rotation flag
     * @param funcPos Function's position in transformation matrices
     * @return        double: function value given 
     */
    private double ellipsFunc (double [] x, int nx, int s_flag, int r_flag, int funcPos) {
        double result = 0.0;
        shiftRotateFunc(x, z, 1.0, s_flag, r_flag, funcPos);
        for (int i = 0; i < nx; i++)
            result += Math.pow(10.0, 6.0 * i / (nx - 1)) * z[i] * z[i];
        return result;
    } // ellipsFunc

    /**
     * Discus Function.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param s_flag  Function shift flag
     * @param r_flag  Function rotation flag
     * @param funcPos Function's position in transformation matrices
     * @return        double: function value given 
     */
    private double discusFunc (double [] x, int nx, int s_flag, int r_flag, int funcPos) {
        shiftRotateFunc(x, z, 1.0, s_flag, r_flag, funcPos); // Shift & rotate
        double result = Math.pow(10.0,6.0) * z[0] * z[0];
        for (int i = 1; i < nx; i++)
            result += z[i] * z[i];
        return result;
    } // discusFunc

    /**
     * Ackley's Function.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param s_flag  Function shift flag
     * @param r_flag  Function rotation flag
     * @param funcPos Function's position in transformation matrices
     * @return        double: function value given 
     */
    private double ackleyFunc (double [] x, int nx, int s_flag, int r_flag, int funcPos) {
        double result = 0.0, sum1 = 0.0, sum2 = 0.0;
        shiftRotateFunc(x, z, 1.0, s_flag, r_flag, funcPos); // Shift & rotate
        for (int i = 0; i < nx; i++) {
            sum1 += z[i] * z[i];
            sum2 += Math.cos(2.0 * Math.PI * z[i]);
        }
        sum1 = -0.2 * Math.sqrt(sum1 / nx);
	    sum2 /= nx;
        result = Math.E - 20.0 * Math.exp(sum1) - Math.exp(sum2) + 20.0;
        return result;
    } // ackleyFunc

    /**
     * Weierstrass Function.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param s_flag  Function shift flag
     * @param r_flag  Function rotation flag
     * @param funcPos Function's position in transformation matrices
     * @return        double: function value given 
     */
    private double weierstrassFunc (double [] x, int nx, int s_flag, int r_flag, int funcPos) {
        int k_max = 20;
        double result = 0.0, sum = 0.0, sum2 = 0.0, a = 0.5, b = 3.0;
        shiftRotateFunc(x, z, 0.5 / 100.0, s_flag, r_flag, funcPos); // Shift & rotate
        for (int i = 0; i < nx; i++) {
            sum = 0.0;
            sum2 = 0.0;
            for (int j = 0; j <= k_max; j++) {
                sum += Math.pow(a, j) * Math.cos(2.0 * Math.PI * Math.pow(b, j) * (z[i] + 0.5));
                sum2 += Math.pow(a, j) * Math.cos(2.0 * Math.PI * Math.pow(b, j) * 0.5);
            }
            result += sum;
        }
        result -= nx * sum2;
        return result;
    } // weierstrassFunc

    /**
     * Griewank's Function.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param s_flag  Function shift flag
     * @param r_flag  Function rotation flag
     * @param funcPos Function's position in transformation matrices
     * @return        double: function value given 
     */
    private double griewankFunc (double [] x, int nx, int s_flag, int r_flag, int funcPos) {
        double result = 0.0, s = 0.0, p = 1.0;
        shiftRotateFunc(x, z, 600.0 / 100.0, s_flag, r_flag, funcPos); // Shift & rotate
        for (int i = 0; i < nx; i++) {
            s += z[i] * z[i];
            p *= Math.cos(z[i] / Math.sqrt(1.0 + i));
        }
        result = 1.0 + s / 4000.0 - p;
        return result;
    } // griewankFunc

    /**
     * Katsuura Function.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param s_flag  Function shift flag
     * @param r_flag  Function rotation flag
     * @param funcPos Function's position in transformation matrices
     * @return        double: function value given 
     */
    private double katsuuraFunc (double [] x, int nx, int s_flag, int r_flag, int funcPos) {
        double temp,tmp1,tmp2,tmp3;
        double result = 1.0;
        tmp3 = Math.pow(1.0 * nx,1.2);
        shiftRotateFunc(x, z, 5.0 / 100.0, s_flag, r_flag, funcPos); // Shift & rotate
        for (int i = 0; i < nx; i++) {
            temp = 0.0;
            for (int j = 1; j <= 32; j++) {
                tmp1 = Math.pow(2.0, j);
                tmp2 = tmp1 * z[i];
                temp += Math.abs(tmp2 - Math.floor(tmp2 + 0.5)) / tmp1;
            }
            result *= Math.pow(1.0 + (i + 1) * temp, 10.0 / tmp3);
        }
        tmp1 = 10.0 / nx / nx;
        result = result * tmp1 - tmp1;
        return result;
    } // katsuuraFunc

    /**
     * Expanded Griewank's plus Rosenbrock's Function.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param s_flag  Function shift flag
     * @param r_flag  Function rotation flag
     * @param funcPos Function's position in transformation matrices
     * @return        double: function value given 
     */
    private double grieRosenFunc (double [] x, int nx, int s_flag, int r_flag, int funcPos) {
        double temp,tmp1,tmp2;
        double result = 0.0;
        shiftRotateFunc(x, z, 5.0 / 100.0, s_flag, r_flag, funcPos); // Shift & rotate
        z[0] += 1.0; // "Shift to origin"
        for (int i = 0; i < nx - 1; i++) {
            z[i + 1] += 1.0; // "Shift to origin"
            tmp1 = z[i] * z[i] - z[i + 1];
            tmp2 = z[i] - 1.0;
            temp = 100.0 * tmp1 * tmp1 + tmp2 * tmp2;
            result += (temp * temp) / 4000.0 - Math.cos(temp) + 1.0; 
        }
        tmp1 = z[nx - 1] * z[nx - 1] - z[0];
        tmp2 = z[nx - 1] - 1.0;
        temp = 100.0 * tmp1 * tmp1 + tmp2 * tmp2;
        result += (temp * temp) / 4000.0 - Math.cos(temp) + 1.0;
        return result;
    } // grieRosenFunc

    /**
     * Expanded Schaffer's F6 Function.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param s_flag  Function shift flag
     * @param r_flag  Function rotation flag
     * @param funcPos Function's position in transformation matrices
     * @return        double: function value given 
     */
    private double eschaffer6Func (double [] x, int nx, int s_flag, int r_flag, int funcPos) {
        double result = 0.0;
        double temp1, temp2;
        shiftRotateFunc(x, z, 1.0, s_flag, r_flag, funcPos); // Shift & rotate
        for (int i = 0; i < nx - 1; i++) {
            temp1 = Math.sin(Math.sqrt(z[i] * z[i] + z[i + 1] * z[i + 1]));
            temp1 = temp1 * temp1;
            temp2 = 1.0 + 0.001 * (z[i] * z[i] + z[i + 1] * z[i + 1]);
            result += 0.5 + (temp1 - 0.5) / (temp2 * temp2);
        }
        temp1 = Math.sin(Math.sqrt(z[nx - 1] * z[nx - 1] + z[0] * z[0]));
        temp1 = temp1 * temp1;
        temp2 = 1.0 + 0.001 * (z[nx - 1] * z[nx - 1] + z[0] * z[0]);
        result += 0.5 + (temp1 - 0.5) / (temp2 * temp2);
        return result;
    } // eschaffer6Func

    /**
     * HappyCat Function.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param s_flag  Function shift flag
     * @param r_flag  Function rotation flag
     * @param funcPos Function's position in transformation matrices
     * @return        double: function value given 
     */
    private double happyCatFunc (double [] x, int nx, int s_flag, int r_flag, int funcPos) {
        // "Original global optimum [-1, -1, ..., -1]"
        double result = 0.0, alpha = 1.0 / 8.0, r2 = 0.0, sum_z = 0.0;
        shiftRotateFunc(x, z, 5.0 / 100.0, s_flag, r_flag, funcPos); // Shift & rotate
        for (int i = 0; i < nx; i++) {
            z[i] = z[i] - 1.0; // "Shift to origin"
            r2 += z[i] * z[i];
            sum_z += z[i];
        }
        result = Math.pow(Math.abs(r2 - nx), 2 * alpha) + (0.5 * r2 + sum_z) / nx + 0.5;
        return result;
    } // happyCatFunc

    /**
     * HGBat Function.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param s_flag  Function shift flag
     * @param r_flag  Function rotation flag
     * @param funcPos Function's position in transformation matrices
     * @return        double: function value given 
     */
    private double hgBatFunc (double [] x, int nx, int s_flag, int r_flag, int funcPos) {
        // "Original global optimum [-1, -1, ..., -1]"
        double result = 0.0, alpha = 1.0 / 4.0, r2 = 0.0, sum_z = 0.0;
        shiftRotateFunc(x, z, 5.0 / 100.0, s_flag, r_flag, funcPos); // Shift & rotate
        for (int i = 0; i < nx; i++) {
            z[i] = z[i] - 1.0; // "Shift to origin"
            r2 += z[i] * z[i];
            sum_z += z[i];
        }
        result = Math.pow(Math.abs(Math.pow(r2, 2.0) - Math.pow(sum_z, 2.0)), 2 * alpha) + (0.5 * r2 + sum_z) / nx + 0.5;
        return result;
    } // hgBatFunc

    /**
     * F1: Shifted and Rotated Bent Cigar Function.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param s_flag  Function shift flag
     * @param r_flag  Function rotation flag
     * @param funcPos Function's position in transformation matrices
     * @return        double: function value given 
     */
    private double bentCigarFunc (double [] x, int nx, int s_flag, int r_flag, int funcPos) {
        double result = 0.0;
        shiftRotateFunc(x, z, 1.0, s_flag, r_flag, funcPos); // Shift & rotate
        result = z[0] * z[0];
        for (int i = 1; i < nx; i++)
            result += Math.pow(10.0, 6.0) * z[i] * z[i];
        return result;
    } // bentCigarFunc

    /**
     * F2: Shifted and Rotated Zakharov Function.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param s_flag  Function shift flag
     * @param r_flag  Function rotation flag
     * @param funcPos Function's position in transformation matrices
     * @return        double: function value given 
     */
    private double zakharovFunc (double [] x, int nx, int s_flag, int r_flag, int funcPos) {
        double result = 0.0;
        double sum1 = 0.0, sum2 = 0.0;
        shiftRotateFunc(x, z, 1.0, s_flag, r_flag, funcPos); // Shift & rotate
        for (int i = 0; i < nx; i++) {
            sum1 += Math.pow(z[i], 2.0);
            sum2 += 0.5 * (i + 1) * z[i];
        }
        result = sum1 + Math.pow(sum2, 2) + Math.pow(sum2, 4);
        return result;
    } // zakharovFunc

    /**
     * F3: Shifted and Rotated Rosenbrock’s Function.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param s_flag  Function shift flag
     * @param r_flag  Function rotation flag
     * @param funcPos Function's position in transformation matrices
     * @return        double: function value given 
     */
    private double rosenbrockFunc (double [] x, int nx, int s_flag, int r_flag, int funcPos) {
        double result = 0.0;
        double temp1, temp2;
        shiftRotateFunc(x, z, 2.048 / 100.0, s_flag, r_flag, funcPos); // Shift & rotate
        z[0] += 1.0; // Shift to origin
        for (int i = 0; i < nx - 1; i++) {
            z[i + 1] += 1.0; // Shift to origin
            temp1 = z[i] * z[i] - z[i + 1];
            temp2 = z[i] - 1.0;
            result += 100.0 * temp1 * temp1 + temp2 * temp2;
        }
        return result;
    } // rosenbrockFunc

    /**
     * F4: Shifted and Rotated Rastrigin’s Function.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param s_flag  Function shift flag
     * @param r_flag  Function rotation flag
     * @param funcPos Function's position in transformation matrices
     * @return        double: function value given 
     */
    private double rastriginFunc (double [] x, int nx, int s_flag, int r_flag, int funcPos) {
        double result = 0.0;
        shiftRotateFunc(x, z, 5.12 / 100.0, s_flag, r_flag, funcPos); // Shift & rotate
        for (int i = 0; i < nx; i++)
            result += (z[i] * z[i] - 10.0 * Math.cos(2.0 * Math.PI * z[i]) + 10.0);
        return result;
    } // rastriginFunc

    /**
     * F5: Shifted and Rotated Schaffer F7 Function - Schwefel's 1.2.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param s_flag  Function shift flag
     * @param r_flag  Function rotation flag
     * @param funcPos Function's position in transformation matrices
     * @return        double: function value given 
     */
    private double schafferF7Func (double [] x, int nx, int s_flag, int r_flag, int funcPos) {
        double result = 0.0;
        double temp;
        shiftRotateFunc(x, z, 1.0, s_flag, r_flag, funcPos); // Shift & rotate
        for (int i = 0; i < nx - 1; i++) {
            z[i] = Math.pow(y[i] * y[i] + y[i + 1] * y[i + 1], 0.5);
            temp = Math.sin(50.0 * Math.pow(z[i], 0.2));
            result += Math.pow(z[i], 0.5) + Math.pow(z[i], 0.5) * temp * temp;
        }
        result = result * result / (nx - 1) / (nx - 1);
        return result;
    } // schafferF7Func

    /**
     * F6: Shifted and Rotated Lunacek Bi-Rastrigin’s Function.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param s_flag  Function shift flag
     * @param r_flag  Function rotation flag
     * @param funcPos Function's position in transformation matrices
     * @return        double: function value given 
     */
    private double biRastriginFunc (double [] x, int nx, int s_flag, int r_flag, int funcPos) {
        double result = 0.0;
        double mu0 = 2.5, d = 1.0, s, mu1, temp, temp1, temp2;
        double [] tempX = new double [nx];
        s = 1.0 - 1.0 / (2.0 * Math.pow(nx + 20.0, 0.5) - 8.2);
	    mu1 = -Math.pow((mu0 * mu0 - d) / s, 0.5);
        if (s_flag == 1)
            shiftFunc(x, y, funcPos);
        else
            for (int i = 0; i < nx; i++) 
                y[i] = x[i]; // "Shrink to original search range"
        for (int i = 0; i < nx; i++) // "Shrink to original search range"
            y[i] *= 10.0 / 100.0;
        for (int i = 0; i < nx; i++) {
            tempX[i] = 2.0 * y[i];
            if (oShift[funcPos][i] < 0.0)
                tempX[i] *= -1.0;
        }
        for (int i = 0; i < nx; i++) {
            z[i] = tempX[i];
            tempX[i] += mu0;
        }
        temp1 = 0.0; temp2 = 0.0;
        for (int i = 0; i < nx; i++) {
            temp = tempX[i] - mu0;
            temp1 += temp * temp;
            temp = tempX[i] - mu1;
            temp2 += temp * temp;
        }
        temp2 *= s;
        temp2 += d * nx;
        temp = 0.0;
        y = new double [nx]; // Need to reset y if it already gets used from a shift
        if (r_flag == 1) {
            rotateFunc(z, y, funcPos);
            for (int i = 0; i < nx; i++)
                temp += Math.cos(2.0 * Math.PI * y[i]);
            if (temp1 < temp2) result = temp1;
            else result = temp2;
            result += 10.0 * (nx - temp);
        }
        else {
            for (int i = 0; i < nx; i++)
                temp += Math.cos(2.0 * Math.PI * z[i]);
            if (temp1 < temp2) result = temp1;
            else result = temp2;
            result += 10.0 * (nx - temp);
        }
        return result;
    } // biRastriginFunc

    /**
     * F7: Shifted and Rotated Non-Continuous Rastrigin’s Function.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param s_flag  Function shift flag
     * @param r_flag  Function rotation flag
     * @param funcPos Function's position in transformation matrices
     * @return        double: function value given 
     */
    private double stepRastriginFunc (double [] x, int nx, int s_flag, int r_flag, int funcPos) {
        double result = 0.0;
        for (int i = 0; i < nx; i++)
            if (Math.abs(y[i] - oShift[0][i]) > 0.5)
                y[i] = oShift[0][i] + Math.floor(2 * (y[i] - oShift[0][i]) + 0.5) / 2.0;
        shiftRotateFunc(x, z, 5.12 / 100.0, s_flag, r_flag, funcPos); // Shift & rotate 
        for (int i = 0; i < nx; i++)
            result += (z[i] * z[i] - 10.0 * Math.cos(2.0 * Math.PI * z[i]) + 10.0);
        return result;
    } // stepRastriginFunc

    /**
     * F8: Shifted and Rotated LevyFunction.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param s_flag  Function shift flag
     * @param r_flag  Function rotation flag
     * @param funcPos Function's position in transformation matrices
     * @return        double: function value given 
     */
    private double levyFunc (double [] x, int nx, int s_flag, int r_flag, int funcPos) {
        double result = 0.0;
        shiftRotateFunc(x, z, 1.0, s_flag, r_flag, funcPos); // Shift & rotate
        double [] w = new double [nx];
        for (int i = 0; i < nx; i++)
            w[i] = 1.0 + (z[i] - 1.0) / 4.0;
        double term1 = Math.pow((Math.sin(Math.PI * w[0])), 2);
        double term3 = Math.pow((w[nx - 1] - 1),2) * (1 + Math.pow((Math.sin(2 * Math.PI* w [nx - 1])),2));
        double sum = 0.0;
        for (int i = 0; i < nx - 1; i++)
            sum += Math.pow((w[i] - 1),2) * (1 + 10 * Math.pow((Math.sin(Math.PI * w[i] + 1)),2));
        result = term1 + sum + term3;
        return result;
    } // levyFunc

    /**
     * F9: Shifted and Rotated Schwefel’s Function.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param s_flag  Function shift flag
     * @param r_flag  Function rotation flag
     * @param funcPos Function's position in transformation matrices
     * @return        double: function value given 
     */
    private double schwefelFunc (double [] x, int nx, int s_flag, int r_flag, int funcPos) {
        double result = 0.0, temp;
        shiftRotateFunc(x, z, 1000.0 / 100.0, s_flag, r_flag, funcPos); // Shift & rotate
        for (int i = 0; i < nx; i++) {
            z[i] += 4.209687462275036e+002;
            if (z[i] > 500) {
                result -= (500.0 - fmod(z[i], 500)) * Math.sin(Math.pow(500.0 - fmod(z[i], 500),0.5));
                temp = (z[i] - 500.0) / 100.0;
                result += temp * temp / nx;
            }
            else if (z[i] < -500) {
                result -= (-500.0 + fmod(Math.abs(z[i]), 500)) * Math.sin(Math.pow(500.0 - fmod(Math.abs(z[i]), 500),0.5));
                temp = (z[i] + 500.0) / 100.0;
                result += temp * temp / nx;
            }
            else
                result -= z[i] * Math.sin(Math.pow(Math.abs(z[i]),0.5));
        }
        result += 4.189828872724338e+002 * nx;
        return result;
    } // schwefelFunc

    /**
    * Computes the floating-point remainder of a/b.
    * @param a Numerator  
    * @param b Denominator
    * @return  double: a % b
    */
    private double fmod(double a, double b) {
        int result = (int) Math.floor(a / b);
        return a - result * b;
    } // fmod

    /**
     * F10: Hybrid function 1.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param s_flag  Function shift flag
     * @param r_flag  Function rotation flag
     * @param funcPos Function's position in transformation matrices
     * @return        double: function value given x
     */
    private double hf01 (double [] x, int nx, int s_flag, int r_flag, int funcPos) {
        double result = 0.0;
        int i, tmp = 0, cf_num = 3;
        int [] G = new int [3];
        int [] G_n = new int [3];
        double [] Gp = {0.2, 0.4, 0.4};

        // Number of indices in each sub group
        for (i = 0; i < cf_num - 1; i++) {
            G_n[i] = (int) Math.ceil(Gp[i] * nx);
            tmp += G_n[i];
        }
        G_n[cf_num - 1] = nx - tmp;
        G[0] = 0;
        for (i = 1; i < cf_num; i++)
            G[i] = G[i - 1] + G_n[i - 1];

        shiftRotateFunc(x, z, 1.0, s_flag, r_flag, funcPos); // Shift & rotate

        // Shuffle objective indices
        for (i = 0; i < nx; i++)
            y[i] = z[shuffle[funcPos][i] - 1];

        // Creating function substitution 2D array
        double [][] subs = new double [G.length][];
        int j;
        for (i = 0; i < subs.length; i++)
            subs[i] = new double [G_n[i]];
        int count = 0;
        for (i = 0; i < subs.length; i++) {
            for (j = 0; j < subs[i].length; j++) {
                subs[i][j] = y[count];
                count++;
            }
        }

        result += zakharovFunc(subs[0], G_n[0], 0,0, 0);
        result += rosenbrockFunc(subs[1], G_n[1], 0,0, 0);
        result += rastriginFunc(subs[2], G_n[2], 0,0, 0);

        return result;
    } // hf01

    /**
     * F11: Hybrid function 2.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param s_flag  Function shift flag
     * @param r_flag  Function rotation flag
     * @param funcPos Function's position in transformation matrices
     * @return        double: function value given x
     */
    private double hf02 (double [] x, int nx, int s_flag, int r_flag, int funcPos) {
        double result = 0.0;
        int i, tmp = 0, cf_num = 3;
        int [] G = new int [3];
        int [] G_n = new int [3];
        double [] Gp = {0.3, 0.3, 0.4}; // Sum to 1

        // Number of indices in each sub group
        for (i = 0; i < cf_num - 1; i++) {
            G_n[i] = (int) Math.ceil(Gp[i] * nx);
            tmp += G_n[i];
        }
        G_n[cf_num - 1] = nx - tmp;
        G[0] = 0;
        for (i = 1; i < cf_num; i++)
            G[i] = G[i - 1] + G_n[i - 1];

        shiftRotateFunc(x, z, 1.0, s_flag, r_flag, funcPos); // Shift & rotate

        // Shuffle objective indices
        for (i = 0; i < nx; i++)
            y[i] = z[shuffle[funcPos][i] - 1];

        // Creating function substitution 2D array
        double [][] subs = new double [G.length][];
        int j;
        for (i = 0; i < subs.length; i++)
            subs[i] = new double [G_n[i]];
        int count = 0;
        for (i = 0; i < subs.length; i++) {
            for (j = 0; j < subs[i].length; j++) {
                subs[i][j] = y[count];
                count++;
            }
        }

        result += ellipsFunc(subs[0], G_n[0], 0,0, 0);
        result += schwefelFunc(subs[1], G_n[1], 0,0, 0);
        result += bentCigarFunc(subs[2], G_n[2], 0,0, 0);

        return result;
    } // hf02

    /**
     * F12: Hybrid function 3.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param s_flag  Function shift flag
     * @param r_flag  Function rotation flag
     * @param funcPos Function's position in transformation matrices
     * @return        double: function value given x
     */
    private double hf03 (double [] x, int nx, int s_flag, int r_flag, int funcPos) {
        double result = 0.0;
        int i, tmp = 0, cf_num = 3;
        int [] G = new int [3];
        int [] G_n = new int [3];
        double [] Gp = {0.3, 0.3, 0.4}; // Sum to 1

        // Number of indices in each sub group
        for (i = 0; i < cf_num - 1; i++) {
            G_n[i] = (int) Math.ceil(Gp[i] * nx);
            tmp += G_n[i];
        }
        G_n[cf_num - 1] = nx - tmp;
        G[0] = 0;
        for (i = 1; i < cf_num; i++)
            G[i] = G[i - 1] + G_n[i - 1];

        shiftRotateFunc(x, z, 1.0, s_flag, r_flag, funcPos); // Shift & rotate

        // Shuffle objective indices
        for (i = 0; i < nx; i++)
            y[i] = z[shuffle[funcPos][i] - 1];

        // Creating function substitution 2D array
        double [][] subs = new double [G.length][];
        int j;
        for (i = 0; i < subs.length; i++)
            subs[i] = new double [G_n[i]];
        int count = 0;
        for (i = 0; i < subs.length; i++) {
            for (j = 0; j < subs[i].length; j++) {
                subs[i][j] = y[count];
                count++;
            }
        }

        result += bentCigarFunc(subs[0], G_n[0], 0,0, 0);
        result += rosenbrockFunc(subs[1], G_n[1], 0,0, 0);
        result += biRastriginFunc(subs[2], G_n[2], 0,0, 0);

        return result;
    } // hf03

    /**
     * F13: Hybrid function 4.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param s_flag  Function shift flag
     * @param r_flag  Function rotation flag
     * @param funcPos Function's position in transformation matrices
     * @return        double: function value given x
     */
    private double hf04 (double [] x, int nx, int s_flag, int r_flag, int funcPos) {
        double result = 0.0;
        int i, tmp = 0, cf_num = 4;
        int [] G = new int [4];
        int [] G_n = new int [4];
        double [] Gp = {0.2, 0.2, 0.2, 0.4}; // Sum to 1

        // Number of indices in each sub group
        for (i = 0; i < cf_num - 1; i++) {
            G_n[i] = (int) Math.ceil(Gp[i] * nx);
            tmp += G_n[i];
        }
        G_n[cf_num - 1] = nx - tmp;
        G[0] = 0;
        for (i = 1; i < cf_num; i++)
            G[i] = G[i - 1] + G_n[i - 1];

        shiftRotateFunc(x, z, 1.0, s_flag, r_flag, funcPos); // Shift & rotate

        // Shuffle objective indices
        for (i = 0; i < nx; i++)
            y[i] = z[shuffle[funcPos][i] - 1];

        // Creating function substitution 2D array
        double [][] subs = new double [G.length][];
        int j;
        for (i = 0; i < subs.length; i++)
            subs[i] = new double [G_n[i]];
        int count = 0;
        for (i = 0; i < subs.length; i++) {
            for (j = 0; j < subs[i].length; j++) {
                subs[i][j] = y[count];
                count++;
            }
        }

        result += ellipsFunc(subs[0], G_n[0], 0,0, 0);
        result += ackleyFunc(subs[1], G_n[1], 0,0, 0);
        result += schafferF7Func(subs[2], G_n[2], 0,0, 0);
        result += rastriginFunc(subs[3], G_n[3], 0,0, 0);

        return result;
    } // hf04

    /**
     * F14: Hybrid function 5.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param s_flag  Function shift flag
     * @param r_flag  Function rotation flag
     * @param funcPos Function's position in transformation matrices
     * @return        double: function value given x
     */
    private double hf05 (double [] x, int nx, int s_flag, int r_flag, int funcPos) {
        double result = 0.0;
        int i, tmp = 0, cf_num = 4;
        int [] G = new int [4];
        int [] G_n = new int [4];
        double [] Gp = {0.2, 0.2, 0.3, 0.3}; // Sum to 1

        // Number of indices in each sub group
        for (i = 0; i < cf_num - 1; i++) {
            G_n[i] = (int) Math.ceil(Gp[i] * nx);
            tmp += G_n[i];
        }
        G_n[cf_num - 1] = nx - tmp;
        G[0] = 0;
        for (i = 1; i < cf_num; i++)
            G[i] = G[i - 1] + G_n[i - 1];

        shiftRotateFunc(x, z, 1.0, s_flag, r_flag, funcPos); // Shift & rotate

        // Shuffle objective indices
        for (i = 0; i < nx; i++)
            y[i] = z[shuffle[funcPos][i] - 1];

        // Creating function substitution 2D array
        double [][] subs = new double [G.length][];
        int j;
        for (i = 0; i < subs.length; i++)
            subs[i] = new double [G_n[i]];
        int count = 0;
        for (i = 0; i < subs.length; i++) {
            for (j = 0; j < subs[i].length; j++) {
                subs[i][j] = y[count];
                count++;
            }
        }

        result += bentCigarFunc(subs[0], G_n[0], 0,0, 0);
        result += hgBatFunc(subs[1], G_n[1], 0,0, 0);
        result += rastriginFunc(subs[2], G_n[2], 0,0, 0);
        result += rosenbrockFunc(subs[3], G_n[3], 0,0, 0);

        return result;
    } // hf05

    /**
     * F15: Hybrid function 6.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param s_flag  Function shift flag
     * @param r_flag  Function rotation flag
     * @param funcPos Function's position in transformation matrices
     * @return        double: function value given x
     */
    private double hf06 (double [] x, int nx, int s_flag, int r_flag, int funcPos) {
        double result = 0.0;
        int i, tmp = 0, cf_num = 4;
        int [] G = new int [4];
        int [] G_n = new int [4];
        double [] Gp = {0.2, 0.2, 0.3, 0.3}; // Sum to 1

        // Number of indices in each sub group
        for (i = 0; i < cf_num - 1; i++) {
            G_n[i] = (int) Math.ceil(Gp[i] * nx);
            tmp += G_n[i];
        }
        G_n[cf_num - 1] = nx - tmp;
        G[0] = 0;
        for (i = 1; i < cf_num; i++)
            G[i] = G[i - 1] + G_n[i - 1];

        shiftRotateFunc(x, z, 1.0, s_flag, r_flag, funcPos); // Shift & rotate

        // Shuffle objective indices
        for (i = 0; i < nx; i++)
            y[i] = z[shuffle[funcPos][i] - 1];

        // Creating function substitution 2D array
        double [][] subs = new double [G.length][];
        int j;
        for (i = 0; i < subs.length; i++)
            subs[i] = new double [G_n[i]];
        int count = 0;
        for (i = 0; i < subs.length; i++) {
            for (j = 0; j < subs[i].length; j++) {
                subs[i][j] = y[count];
                count++;
            }
        }

        result += eschaffer6Func(subs[0], G_n[0], 0,0, 0);
        result += hgBatFunc(subs[1], G_n[1], 0,0, 0);
        result += rosenbrockFunc(subs[2], G_n[2], 0,0, 0);
        result += schwefelFunc(subs[3], G_n[3], 0,0, 0);

        return result;
    } // hf06

    /**
     * F16: Hybrid function 7.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param s_flag  Function shift flag
     * @param r_flag  Function rotation flag
     * @param funcPos Function's position in transformation matrices
     * @return        double: function value given x
     */
    private double hf07 (double [] x, int nx, int s_flag, int r_flag, int funcPos) {
        double result = 0.0;
        int i, tmp = 0, cf_num = 5;
        int [] G = new int [5];
        int [] G_n = new int [5];
        double [] Gp = {0.1, 0.2, 0.2, 0.2, 0.3}; // Sum to 1

        // Number of indices in each sub group
        for (i = 0; i < cf_num - 1; i++) {
            G_n[i] = (int) Math.ceil(Gp[i] * nx);
            tmp += G_n[i];
        }
        G_n[cf_num - 1] = nx - tmp;
        G[0] = 0;
        for (i = 1; i < cf_num; i++)
            G[i] = G[i - 1] + G_n[i - 1];

        shiftRotateFunc(x, z, 1.0, s_flag, r_flag, funcPos); // Shift & rotate

        // Shuffle objective indices
        for (i = 0; i < nx; i++)
            y[i] = z[shuffle[funcPos][i] - 1];

        // Creating function substitution 2D array
        double [][] subs = new double [G.length][];
        int j;
        for (i = 0; i < subs.length; i++)
            subs[i] = new double [G_n[i]];
        int count = 0;
        for (i = 0; i < subs.length; i++) {
            for (j = 0; j < subs[i].length; j++) {
                subs[i][j] = y[count];
                count++;
            }
        }

        result += katsuuraFunc(subs[0], G_n[0], 0,0, 0);
        result += ackleyFunc(subs[1], G_n[1], 0,0, 0);
        result += grieRosenFunc(subs[2], G_n[2], 0,0, 0);
        result += schwefelFunc(subs[3], G_n[3], 0,0, 0);
        result += rastriginFunc(subs[4], G_n[4], 0,0, 0);

        return result;
    } // hf07

    /**
     * F17: Hybrid function 8.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param s_flag  Function shift flag
     * @param r_flag  Function rotation flag
     * @param funcPos Function's position in transformation matrices
     * @return        double: function value given x
     */
    private double hf08 (double [] x, int nx, int s_flag, int r_flag, int funcPos) {
        double result = 0.0;
        int i, tmp = 0, cf_num = 5;
        int [] G = new int [5];
        int [] G_n = new int [5];
        double [] Gp = {0.2, 0.2, 0.2, 0.2, 0.2}; // Sum to 1

        // Number of indices in each sub group
        for (i = 0; i < cf_num - 1; i++) {
            G_n[i] = (int) Math.ceil(Gp[i] * nx);
            tmp += G_n[i];
        }
        G_n[cf_num - 1] = nx - tmp;
        G[0] = 0;
        for (i = 1; i < cf_num; i++)
            G[i] = G[i - 1] + G_n[i - 1];

        shiftRotateFunc(x, z, 1.0, s_flag, r_flag, funcPos); // Shift & rotate

        // Shuffle objective indices
        for (i = 0; i < nx; i++)
            y[i] = z[shuffle[funcPos][i] - 1];

        // Creating function substitution 2D array
        double [][] subs = new double [G.length][];
        int j;
        for (i = 0; i < subs.length; i++)
            subs[i] = new double [G_n[i]];
        int count = 0;
        for (i = 0; i < subs.length; i++) {
            for (j = 0; j < subs[i].length; j++) {
                subs[i][j] = y[count];
                count++;
            }
        }

        result += ellipsFunc(subs[0], G_n[0], 0,0, 0);
        result += ackleyFunc(subs[1], G_n[1], 0,0, 0);
        result += rastriginFunc(subs[2], G_n[2], 0,0, 0);
        result += hgBatFunc(subs[3], G_n[3], 0,0, 0);
        result += discusFunc(subs[4], G_n[4], 0,0, 0);

        return result;
    } // hf08

    /**
     * F18: Hybrid function 9.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param s_flag  Function shift flag
     * @param r_flag  Function rotation flag
     * @param funcPos Function's position in transformation matrices
     * @return        double: function value given x
     */
    private double hf09 (double [] x, int nx, int s_flag, int r_flag, int funcPos) {
        double result = 0.0;
        int i, tmp = 0, cf_num = 5;
        int [] G = new int [5];
        int [] G_n = new int [5];
        double [] Gp = {0.2, 0.2, 0.2, 0.2, 0.2}; // Sum to 1

        // Number of indices in each sub group
        for (i = 0; i < cf_num - 1; i++) {
            G_n[i] = (int) Math.ceil(Gp[i] * nx);
            tmp += G_n[i];
        }
        G_n[cf_num - 1] = nx - tmp;
        G[0] = 0;
        for (i = 1; i < cf_num; i++)
            G[i] = G[i - 1] + G_n[i - 1];

        shiftRotateFunc(x, z, 1.0, s_flag, r_flag, funcPos); // Shift & rotate

        // Shuffle objective indices
        for (i = 0; i < nx; i++)
            y[i] = z[shuffle[funcPos][i] - 1];

        // Creating function substitution 2D array
        double [][] subs = new double [G.length][];
        int j;
        for (i = 0; i < subs.length; i++)
            subs[i] = new double [G_n[i]];
        int count = 0;
        for (i = 0; i < subs.length; i++) {
            for (j = 0; j < subs[i].length; j++) {
                subs[i][j] = y[count];
                count++;
            }
        }

        result += bentCigarFunc(subs[0], G_n[0], 0,0, 0);
        result += rastriginFunc(subs[1], G_n[1], 0,0, 0);
        result += grieRosenFunc(subs[2], G_n[2], 0,0, 0);
        result += weierstrassFunc(subs[3], G_n[3], 0,0, 0);
        result += eschaffer6Func(subs[4], G_n[4], 0,0, 0);

        return result;
    } // hf09

    /**
     * F19: Hybrid function 10.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param s_flag  Function shift flag
     * @param r_flag  Function rotation flag
     * @param funcPos Function's position in transformation matrices
     * @return  double: function value given x
     */
    private double hf10 (double [] x, int nx, int s_flag, int r_flag, int funcPos) {
        double result = 0.0;
        int i, tmp = 0, cf_num = 6;
        int [] G = new int [6];
        int [] G_n = new int [6];
        double [] Gp = {0.1, 0.1, 0.2, 0.2, 0.2, 0.2}; // Sum to 1

        // Number of indices in each sub group
        for (i = 0; i < cf_num - 1; i++) {
            G_n[i] = (int) Math.ceil(Gp[i] * nx);
            tmp += G_n[i];
        }
        G_n[cf_num - 1] = nx - tmp;
        G[0] = 0;
        for (i = 1; i < cf_num; i++)
            G[i] = G[i - 1] + G_n[i - 1];

        shiftRotateFunc(x, z, 1.0, s_flag, r_flag, funcPos); // Shift & rotate

        // Shuffle objective indices
        for (i = 0; i < nx; i++)
            y[i] = z[shuffle[funcPos][i] - 1];

        // Creating function substitution 2D array
        double [][] subs = new double [G.length][];
        int j;
        for (i = 0; i < subs.length; i++)
            subs[i] = new double [G_n[i]];
        int count = 0;
        for (i = 0; i < subs.length; i++) {
            for (j = 0; j < subs[i].length; j++) {
                subs[i][j] = y[count];
                count++;
            }
        }

        result += hgBatFunc(subs[0], G_n[0], 0,0, 0);
        result += katsuuraFunc(subs[1], G_n[1], 0,0, 0);
        result += ackleyFunc(subs[2], G_n[2], 0,0, 0);
        result += rastriginFunc(subs[3], G_n[3], 0,0, 0);
        result += schwefelFunc(subs[4], G_n[4], 0,0, 0);
        result += schafferF7Func(subs[5], G_n[5], 0,0, 0);

        return result;
    } // hf10

    /**
     * F20: Composition function 1.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param r_flag  Function rotation flag
     * @return        double: function value given x
     */
    private double cf01 (double [] x, int nx, int r_flag) {
        int i, cf_num = 3;
        double [] fit = new double [cf_num];
        double [] delta = {10, 20, 30};
        double [] bias = {0, 100, 200};
        
        // Need to specify oShift and rotation data per function
        i = 0;
        fit[i] = rosenbrockFunc(x, nx, 1, r_flag, i);
        i = 1;
        fit[i] = ellipsFunc(x, nx, 1, r_flag, i);
        fit[i] = 10000 * fit[i] / 1e+10;
        i = 2;
        fit[i] = rastriginFunc(x, nx, 1, r_flag, i); 

        return cf_cal(x, nx, delta, bias, fit, cf_num);
    } // cf01

    /**
     * F21: Composition function 2.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param r_flag  Function rotation flag
     * @return        double: function value given x
     */
    private double cf02 (double [] x, int nx, int r_flag) {
        int i, cf_num = 3;
        double [] fit = new double [cf_num];
        double [] delta = {10, 20, 30};
        double [] bias = {0, 100, 200};
        
        // Need to specify oShift and rotation data per function
        i = 0;
        fit[i] = rastriginFunc(x, nx, 1, r_flag, i);
        i = 1;
        fit[i] = griewankFunc(x, nx, 1, r_flag, i);
        fit[i] = 1000 * fit[i] / 100.0;
        i = 2;
        fit[i] = schwefelFunc(x, nx, 1, r_flag, i); 

        return cf_cal(x, nx, delta, bias, fit, cf_num);
    } // cf02

    /**
     * F22: Composition function 3.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param r_flag  Function rotation flag
     * @return        double: function value given x
     */
    private double cf03 (double [] x, int nx, int r_flag) {
        int i, cf_num = 4;
        double [] fit = new double [cf_num];
        double [] delta = {10, 20, 30, 40};
        double [] bias = {0, 100, 200, 300};
        
        // Need to specify oShift and rotation data per function
        i = 0;
        fit[i] = rosenbrockFunc(x, nx, 1, r_flag, i);
        i = 1;
        fit[i] = ackleyFunc(x, nx, 1, r_flag, i);
        fit[i] = 1000 * fit[i] / 100.0;
        i = 2;
        fit[i] = schwefelFunc(x, nx, 1, r_flag, i);
        i = 3;
        fit[i] = rastriginFunc(x, nx, 1, r_flag, i);

        return cf_cal(x, nx, delta, bias, fit, cf_num);
    } // cf03

    /**
     * F23: Composition function 4.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param r_flag  Function rotation flag
     * @return        double: function value given x
     */
    private double cf04 (double [] x, int nx, int r_flag) {
        int i, cf_num = 4;
        double [] fit = new double [cf_num];
        double [] delta = {10, 20, 30, 40};
        double [] bias = {0, 100, 200, 300};
        
        // Need to specify oShift and rotation data per function
        i = 0;
        fit[i] = ackleyFunc(x, nx, 1, r_flag, i);
        fit[i] = 1000 * fit[i] / 100.0;
        i = 1;
        fit[i] = ellipsFunc(x, nx, 1, r_flag, i);
        fit[i] = 10000 * fit[i] / 1e+10;
        i = 2;
        fit[i] = griewankFunc(x, nx, 1, r_flag, i);
        fit[i] = 1000 * fit[i] / 100.0;
        i = 3;
        fit[i] = rastriginFunc(x, nx, 1, r_flag, i);

        return cf_cal(x, nx, delta, bias, fit, cf_num);
    } // cf04

    /**
     * F24: Composition function 5.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param r_flag  Function rotation flag
     * @return        double: function value given x
     */
    private double cf05 (double [] x, int nx, int r_flag) {
        int i, cf_num = 5;
        double [] fit = new double [cf_num];
        double [] delta = {10, 20, 30, 40, 50};
        double [] bias = {0, 100, 200, 300, 400};
        
        // Need to specify oShift and rotation data per function
        i = 0;
        fit[i] = rastriginFunc(x, nx, 1, r_flag, i);
        fit[i] = 10000 * fit[i] / 1e+3;
        i = 1;
        fit[i] = happyCatFunc(x, nx, 1, r_flag, i);
        fit[i] = 1000 * fit[i] / 1e+3;
        i = 2;
        fit[i] = ackleyFunc(x, nx, 1, r_flag, i);
        fit[i] = 1000 * fit[i] / 100.0;
        i = 3;
        fit[i] = discusFunc(x, nx, 1, r_flag, i);
        fit[i] = 10000 * fit[i] / 1e+10;	
        i = 4;
        fit[i] = rosenbrockFunc(x, nx, 1, r_flag, i);

        return cf_cal(x, nx, delta, bias, fit, cf_num);
    } // cf05

    /**
     * F25: Composition function 6.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param r_flag  Function rotation flag
     * @return        double: function value given x
     */
    private double cf06 (double [] x, int nx, int r_flag) {
        int i, cf_num = 5;
        double [] fit = new double [cf_num];
        double [] delta = {10, 20, 20, 30, 40};
        double [] bias = {0, 100, 200, 300, 400};
        
        // Need to specify oShift and rotation data per function
        i = 0;
        fit[i] = eschaffer6Func(x, nx, 1, r_flag, i);
        fit[i] = 10000 * fit[i] / 2e+7;
        i = 1;
        fit[i] = schwefelFunc(x, nx, 1, r_flag, i);
        i = 2;
        fit[i] = griewankFunc(x, nx, 1, r_flag, i);
        fit[i] = 1000 * fit[i] / 100.0;
        i = 3;
        fit[i] = rosenbrockFunc(x, nx, 1, r_flag, i);
        i = 4;
        fit[i] = rastriginFunc(x, nx, 1, r_flag, i);
        fit[i] = 10000 * fit[i] / 1e+3;

        return cf_cal(x, nx, delta, bias, fit, cf_num);
    } // cf06

    /**
     * F26: Composition function 7.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param r_flag  Function rotation flag
     * @return        double: function value given x
     */
    private double cf07 (double [] x, int nx, int r_flag) {
        int i, cf_num = 6;
        double [] fit = new double [cf_num];
        double [] delta = {10, 20, 30, 40, 50, 60};
        double [] bias = {0, 100, 200, 300, 400, 500};

        i = 0;
        fit[i] = hgBatFunc(x, nx, 1, r_flag, i);
        fit[i] = 10000 * fit[i] / 1000;
        i = 1;
        fit[i] = rastriginFunc(x, nx, 1, r_flag, i);
        fit[i] = 10000 * fit[i] / 1e+3;
        i = 2;
        fit[i] = schwefelFunc(x, nx, 1, r_flag, i);
        fit[i] = 10000 * fit[i] / 4e+3;
        i = 3;
        fit[i] = bentCigarFunc(x, nx, 1, r_flag, i);
        fit[i] = 10000 * fit[i] / 1e+30;
        i = 4;
        fit[i] = ellipsFunc(x, nx, 1, r_flag, i);
        fit[i] = 10000 * fit[i] / 1e+10;
        i = 5;
        fit[i] = eschaffer6Func(x, nx, 1, r_flag, i);
        fit[i] = 10000 * fit[i] / 2e+7;
 
        return cf_cal(x, nx, delta, bias, fit, cf_num);
    } // cf07

    /**
     * F27: Composition function 8.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param r_flag  Function rotation flag
     * @return        double: function value given x
     */
    private double cf08 (double [] x, int nx, int r_flag) {
        int i, cf_num = 6;
        double [] fit = new double [cf_num];
        double [] delta = {10,20,30,40,50,60};
        double [] bias = {0, 100, 200, 300, 400, 500};
        
        i = 0;
        fit[i] = ackleyFunc(x, nx, 1, r_flag, i);
        fit[i] = 1000 * fit[i] / 100;
        i = 1;
        fit[i] = griewankFunc(x, nx, 1, r_flag, i);
        fit[i] = 1000 * fit[i] / 100;
        i = 2;
        fit[i] = discusFunc(x, nx, 1, r_flag, i);
        fit[i] = 10000 * fit[i] / 1e+10;
        i = 3;
        fit[i] = rosenbrockFunc(x, nx, 1, r_flag, i);
        i = 4;
        fit[i] = happyCatFunc(x, nx, 1, r_flag, i);
        fit[i] = 1000 * fit[i] / 1e+3;
        i = 5;
        fit[i] = eschaffer6Func(x, nx, 1, r_flag, i);
        fit[i] = 10000 * fit[i] / 2e+7;

        return cf_cal(x, nx, delta, bias, fit, cf_num);
    } // cf08

    /**
     * F28: Composition function 9.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param r_flag  Function rotation flag
     * @return        double: function value given x
     */
    private double cf09 (double [] x, int nx, int r_flag) {
        int i, cf_num = 3;
        double [] fit = new double [cf_num];
        double [] delta = {10, 30, 50};
        double [] bias = {0, 100, 200};

        i = 0;
        fit[i] = hf05(x, nx, 1, r_flag, i);
        i = 1;
        fit[i] = hf06(x, nx, 1, r_flag, i);
        i = 2;
        fit[i] = hf07(x, nx, 1, r_flag, i);
        
        return cf_cal(x, nx, delta, bias, fit, cf_num);
    } // cf09

    /**
     * F29: Composition function 10.
     * @param x       Decision vector
     * @param nx      Number of dimensions
     * @param r_flag  Function rotation flag
     * @return        double: function value given x
     */
    private double cf10 (double [] x, int nx, int r_flag) {
        int i, cf_num = 3;
        double [] fit = new double [cf_num];
        double [] delta = {10, 30, 50};
        double [] bias = {0, 100, 200};

        i = 0;
        fit[i] = hf05(x, nx, 1, r_flag, i);
        i = 1;
        fit[i] = hf08(x, nx, 1, r_flag, i);
        i = 2;
        fit[i] = hf09(x, nx, 1, r_flag, i);
        
        return cf_cal(x, nx, delta, bias, fit, cf_num);
    } // cf10

    /**
     * Computes the objective value of the decision vector for this objects function.
     * @param x Decision vector
     * @return  double: objective fitness
     */
    public double compute (double [] x) {

        double f = Double.MAX_VALUE; // Final fitness

        // Temporary vectors - need to be reset each time a decision vector is evaluated
        this.y = new double [x.length];
        this.z = new double [x.length];

        switch (funcNum) {
            case 1: // Shifted and Rotated Bent Cigar Function
                f = bentCigarFunc(x, n, 1, 1, 0);
                return f + 100.0;
            case 2: // Shifted and Rotated Zakharov Function
                f = zakharovFunc(x, n, 1, 1, 0);  
                return f + 200.0;
            case 3: // Shifted and Rotated Rosenbrock’s Function
                f = rosenbrockFunc(x, n, 1, 1, 0); 
                return f + 300.0;
            case 4: // Shifted and Rotated Rastrigin’s Function
                f = rastriginFunc(x, n, 1, 1, 0);  
                return f + 400.0;
            case 5: // Shifted and Rotated Schaffer F7 Function
                f = schafferF7Func(x, n, 1, 1, 0); 
                return f + 500.0;
            case 6: // Shifted and Rotated Lunacek Bi-Rastrigin’s Function
                f = biRastriginFunc(x, n, 1, 1, 0); 
                return f + 600.0;
            case 7: // Shifted and Rotated Non-Continuous Rastrigin’s Function
                f = stepRastriginFunc(x, n, 1, 1, 0); 
                return f + 700.0;
            case 8: // Shifted and RotatedLevyFunction
                f = levyFunc(x, n, 1, 1, 0); 
                return f + 800.0;
            case 9: // Shifted and Rotated Schwefel’s Function
                f = schwefelFunc(x, n, 1, 1, 0); 
                return f + 900.0;
            case 10: // Hybrid function 1
                f = hf01(x, n, 1, 1, 0); 
                return f + 1000.0;
            case 11: // Hybrid function 2
                f = hf02(x, n, 1, 1, 0); 
                return f + 1100.0;
            case 12: // Hybrid function 3
                f = hf03(x, n, 1, 1, 0); 
                return f + 1200.0;
            case 13: // Hybrid function 4
                f = hf04(x, n, 1, 1, 0); 
                return f + 1300.0;
            case 14: // Hybrid function 5
                f = hf05(x, n, 1, 1, 0); 
                return f + 1400.0;
            case 15: // Hybrid function 6
                f = hf06(x, n, 1, 1, 0); 
                return f + 1500.0;
            case 16: // Hybrid function 7
                f = hf07(x, n, 1, 1, 0); 
                return f + 1600.0;
            case 17: // Hybrid function 8
                f = hf08(x, n, 1, 1, 0); 
                return f + 1700.0;
            case 18: // Hybrid function 9
                f = hf09(x, n, 1, 1, 0); 
                return f + 1800.0;
            case 19: // Hybrid function 10
                f = hf10(x, n, 1, 1, 0); 
                return f + 1900.0;
            case 20: // Composition function 1
                f = cf01(x, n, 1); 
                return f + 2000.0;
            case 21: // Composition function 2
                f = cf02(x, n, 1); 
                return f + 2100.0;
            case 22: // Composition function 3
                f = cf03(x, n, 1); 
                return f + 2200.0;
            case 23: // Composition function 4
                f = cf04(x, n, 1); 
                return f + 2300.0;
            case 24: // Composition function 5
                f = cf05(x, n, 1); 
                return f + 2400.0;
            case 25: // Composition function 6
                f = cf06(x, n, 1); 
                return f + 2500.0;
            case 26: // Composition function 7
                f = cf07(x, n, 1); 
                return f + 2600.0;
            case 27: // Composition function 8
                f = cf08(x, n, 1); 
                return f + 2700.0;
            case 28: // Composition function 9
                f = cf09(x, n, 1); 
                return f + 2800.0;
            case 29: // Composition function 10
                f = cf10(x, n, 1); 
                return f + 2900.0;
            default: // Worst possible case
                return f;
        }

    } // compute

    /**
     * @return double: the global optimum of the function.
     */
    public double getOptimum () {
        return funcNum * 100.0;
    } // getOptimum

    /**
     * @return int: the number associated with this objects objective function.
     */
    public int getFuncNum () {
        return funcNum;
    } // getFuncNum

    /**
     * @return String: the full name of this objects function.
     */
    public String getName () {
        switch (funcNum) {
            case 1:  return "Shifted and Rotated Bent Cigar Function";
            case 2:  return "Shifted and Rotated Zakharov Function";
            case 3:  return "Shifted and Rotated Rosenbrock's Function";
            case 4:  return "Shifted and Rotated Rastrigin's Function";
            case 5:  return "Shifted and Rotated Schaffer F7 Function";
            case 6:  return "Shifted and Rotated Lunacek Bi-Rastrigin's Function";
            case 7:  return "Shifted and Rotated Non-Continuous Rastrigin's Function";
            case 8:  return "Shifted and Rotated Levy Function";
            case 9:  return "Shifted and Rotated Schwefel's Function";
            case 10: return "Hybrid function 1";
            case 11: return "Hybrid function 2";
            case 12: return "Hybrid function 3";
            case 13: return "Hybrid function 4";
            case 14: return "Hybrid function 5";
            case 15: return "Hybrid function 6";
            case 16: return "Hybrid function 7";
            case 17: return "Hybrid function 8";
            case 18: return "Hybrid function 9";
            case 19: return "Hybrid function 10";
            case 20: return "Composition function 1";
            case 21: return "Composition function 2";
            case 22: return "Composition function 3";
            case 23: return "Composition function 4";
            case 24: return "Composition function 5";
            case 25: return "Composition function 6";
            case 26: return "Composition function 7";
            case 27: return "Composition function 8";
            case 28: return "Composition function 9";
            case 29: return "Composition function 10";
            default: return "Function not found";
        }
    } // getName

    /**
     * @return String: the function number with an F in front.
     */
    public String getShortName () {
        return "F" + funcNum;
    } // getShortName

    /**
     * @return double: the minimum bound of the search space.
     */
    public double getMin () {
        return min;
    } // getMin

    /**
     * @return double: the maximum bound of the search space.
     */
    public double getMax () {
        return max;
    } // getMax

    /**
     * @return int: the dimension of the search space.
     */
    public int getDimension () {
        return n;
    } // getDimension

    /**
     * @return double: the absolute difference between the bounds of the search space.
     */
    public double getRange () {
        return range;
    } // getRange

} // Function