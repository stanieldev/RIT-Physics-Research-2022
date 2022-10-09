//
//  Drop.cpp
//  DropletShape
//
//  Created by Kara Maki on 9/17/21.
//  Copyright © 2021 Kara Maki. All rights reserved.
//

#include "Redo.h"
#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#define M_PI 3.14159265358979323846



// store the values, coordinates of points into the arrays.
void Drop::CalculateTimeArrays()
{
    std::cout << "CalculateTimeArrays ran\n";
    int m_nPoints;
    int m_mPoints;
    int m_nTimes;
    m_nPoints = 41; // odd number greater than 7
    m_mPoints = 41; // odd number greater than 7
    m_nTimes = 2101;
    for (int nCount = 0; nCount < m_nTimes; nCount++) //time steps
    {
        for (int i = 0; i < m_nPoints; i++) // i array
        {
            for (int j = 0; j < m_mPoints; j++) // j array
            {
                m_pdXt[i][j][nCount] = 0;
                m_pdYt[i][j][nCount] = 0;
                m_pdZt[i][j][nCount] = 0;
                theta[i][j] = 0;
            }
        }
        std::cout << nCount << "\n";
        double angle; // all needed to be double to stop errors
        double dangle;
        double R;
        double phi;
        double r;
        int i;
        int j;
        double dR;
        angle = 0;
        dangle = (PI / 2) / (m_nPoints - 1); // index from (m_mPoints - 1)
        phi = PI / 4;
        R = 1.9407025E-3;//1E-3;//1.9407025E-1;//1E-3;
        dR = R / (0.5 * (m_mPoints - 1));   // index from (m_mPoints - 1)
        double X;
        double Y;
        //std::cout << m_nTimes <<"\n";
        //std::cout << nCount <<"\n";
        if (nCount == 0)
        {
            std::cout << "if nCount == 0\n";
            std::cout << "boundaries\n";
            // bottom
            std::cout << "bottom\n";
            int j = 0; // index from 1
            for (int i = 0; i < m_nPoints; i++) // index from (i = 1:m_nPoints)
            {
                X = R * cos(angle);
                Y = R * sin(angle);
                m_pdXt[i][j][nCount] = X;
                m_pdYt[i][j][nCount] = Y;
                r = sqrt(X * X + Y * Y);
                m_pdZt[i][j][nCount] = sqrt(((R / sin(phi)) * (R / sin(phi))) - (r * r)) - R / tan(phi);
                if (i == m_nPoints - 1) // index from m_nPoints
                {
                }
                else
                {
                    angle = angle + dangle;
                }
            }
            // right
            std::cout << "right\n";
            i = m_nPoints - 1; // index from m_nPoints
            for (int j = 0; j < m_mPoints; j++) // index from (j = 1:m_nPoints)
            {
                X = R * cos(angle);
                Y = R * sin(angle);
                m_pdXt[i][j][nCount] = X;
                m_pdYt[i][j][nCount] = Y;
                r = sqrt(X * X + Y * Y);
                m_pdZt[i][j][nCount] = sqrt(((R / sin(phi)) * (R / sin(phi))) - (r * r)) - R / tan(phi);
                if (j == m_mPoints - 1) // index from m_nPoints
                {
                }
                else
                {
                    angle = angle + dangle;
                }
            }
            // top
            std::cout << "top\n";
            j = m_nPoints - 1;  // index from m_nPoints
            for (int i = m_nPoints - 1; i > -1; --i) // index from (i = m_nPoints:-1:1);  i=m_nPoints results in an error
            {
                X = R * cos(angle);
                Y = R * sin(angle);
                m_pdXt[i][j][nCount] = X;
                m_pdYt[i][j][nCount] = Y;
                r = sqrt(X * X + Y * Y);
                m_pdZt[i][j][nCount] = sqrt(((R / sin(phi)) * (R / sin(phi))) - (r * r)) - R / tan(phi);
                if (i == 0) // index from 1
                {
                }
                else
                {
                    angle = angle + dangle;
                }
            }
            // left
            std::cout << "left\n";
            i = 0;  // index from 1
            for (int j = m_nPoints - 1; j > -1; --j) // index from (j = m_nPoints:-1:1)        j=m_nPoints results in an error
            {
                X = R * cos(angle);
                Y = R * sin(angle);
                m_pdXt[i][j][nCount] = X;
                m_pdYt[i][j][nCount] = Y;
                r = sqrt(X * X + Y * Y);
                m_pdZt[i][j][nCount] = sqrt(((R / sin(phi)) * (R / sin(phi))) - (r * r)) - R / tan(phi);
                if (j == 0) // index from 1
                {
                }
                else
                {
                    angle = angle + dangle;
                }
            }

            for (int k = 1; k < 0.5 * (m_mPoints - 1); k++) // index from k = 2:0.5*(m_mpoints-1)
            {
                std::cout << "interior\n";
                angle = 0;
                std::cout << k << "\n";
                dangle = (PI / 2) / ((m_nPoints - 1) - 2 * (k)); // index from dangle = (pi/2)/((m_npoints-1)-2*(k-1));
                //std::cout << dangle << "\n";

                // bottom
                std::cout << "interior bottom\n";
                j = k;
                for (int i = k; i < m_nPoints - (k - 1) - 1; i++) // index from k:m_npoints-(k-1)
                {
                    X = (R - (k)*dR) * cos(angle);
                    Y = (R - (k)*dR) * sin(angle);
                    m_pdXt[i][j][nCount] = X;
                    m_pdYt[i][j][nCount] = Y;
                    r = sqrt(X * X + Y * Y);
                    m_pdZt[i][j][nCount] = sqrt(((R / sin(phi)) * (R / sin(phi))) - (r * r)) - R / tan(phi);
                    if (i == (m_nPoints - (k - 1)) - 2) // index from m_npoints-(k-1)
                    {
                    }
                    else
                    {
                        angle = angle + dangle;
                    }
                }

                // right
                std::cout << "interior right\n";
                i = m_nPoints - (k - 1) - 2;  // changed from i=m_nPoints-(k-1) to avoid error
                for (j = k; j < m_nPoints - (k - 1) - 1; j++) // index from k:m_mpoints-(k-1)
                {
                    //std::cout << j << "\n";
                    X = (R - (k)*dR) * cos(angle);
                    Y = (R - (k)*dR) * sin(angle);
                    m_pdXt[i][j][nCount] = X;
                    m_pdYt[i][j][nCount] = Y;
                    r = sqrt(X * X + Y * Y);
                    m_pdZt[i][j][nCount] = sqrt(((R / sin(phi)) * (R / sin(phi))) - (r * r)) - R / tan(phi);
                    if (j == (m_nPoints - (k - 1)) - 2) // index from m_mpoints-(k-1)
                    {
                    }
                    else
                    {
                        angle = angle + dangle;
                    }
                }


                // top
                std::cout << "interior top\n";
                j = m_nPoints - (k - 1) - 2; // changed from j = m_nPoints - (k-1) to avoid errors
                for (int i = m_nPoints - (k - 1) - 2; i > k - 1; i--) // index from m_npoints-(k-1):-1:k       changed to i = m_nPoints - (k - 1) - 1 to avoid error
                {
                    X = (R - (k)*dR) * cos(angle);
                    Y = (R - (k)*dR) * sin(angle);
                    m_pdXt[i][j][nCount] = X;
                    m_pdYt[i][j][nCount] = Y;
                    r = sqrt(X * X + Y * Y);
                    m_pdZt[i][j][nCount] = sqrt(((R / sin(phi)) * (R / sin(phi))) - (r * r)) - R / tan(phi);
                    if (i == k)
                    {
                    }
                    else
                    {
                        angle = angle + dangle;
                    }
                }


                // left
                std::cout << "interior left\n";
                i = k;
                for (j = m_nPoints - (k - 1) - 2; j > k - 1; j--) // index from m_npoints-(k-1):-1:k       changed to j = m_nPoints - (k - 1) - 1 to avoid error
                {
                    //std::cout << j << "\n";
                    X = (R - (k)*dR) * cos(angle);
                    Y = (R - (k)*dR) * sin(angle);
                    m_pdXt[i][j][nCount] = X;
                    m_pdYt[i][j][nCount] = Y;
                    r = sqrt(X * X + Y * Y);
                    m_pdZt[i][j][nCount] = sqrt(((R / sin(phi)) * (R / sin(phi))) - (r * r)) - R / tan(phi);
                    if (j == k)
                    {
                    }
                    else
                    {
                        angle = angle + dangle;
                    }
                }
            }
            std::cout << "set middle\n";
            r = 0;
            m_pdZt[(m_nPoints - 1) / 2][(m_mPoints - 1) / 2][nCount] = sqrt(((R / sin(phi)) * (R / sin(phi))) - r * r) - R / tan(phi);
        }

        else
        {
            for (int i = 0; i < m_nPoints; i++)
            {
                for (int j = 0; j < m_mPoints; j++) {
                    //std::cout << nCount << "nCount \n";
            // bottom left point
                    if (i == 0 && j == 0)
                    {

                        m_pdXt[i][j][nCount] = (m_pdXt[0][1][nCount - 1] + m_pdXt[1][0][nCount - 1]) / 2;
                        m_pdYt[i][j][nCount] = (m_pdYt[0][1][nCount - 1] + m_pdYt[1][0][nCount - 1]) / 2;
                        m_pdZt[i][j][nCount] = 0;
                    }

                    // top left point
                    else if (i == 0 && j == m1)
                    {
                        m_pdXt[i][j][nCount] = (m_pdXt[0][m2][nCount - 1] + m_pdXt[1][m1][nCount - 1]) / 2;
                        m_pdYt[i][j][nCount] = (m_pdYt[0][m2][nCount - 1] + m_pdYt[1][m1][nCount - 1]) / 2;
                        m_pdZt[i][j][nCount] = 0;
                    }

                    // top right point
                    else if (i == n1 && j == m1)
                    {
                        m_pdXt[i][j][nCount] = (m_pdXt[n1][m2][nCount - 1] + m_pdXt[n2][m1][nCount - 1]) / 2;
                        m_pdYt[i][j][nCount] = (m_pdYt[n1][m2][nCount - 1] + m_pdYt[n2][m1][nCount - 1]) / 2;
                        m_pdZt[i][j][nCount] = 0;
                    }

                    // bottom right point
                    else if (i == n1 && j == 0)
                    {
                        m_pdXt[i][j][nCount] = (m_pdXt[n1][1][nCount - 1] + m_pdXt[n2][0][nCount - 1]) / 2;
                        m_pdYt[i][j][nCount] = (m_pdYt[n1][1][nCount - 1] + m_pdYt[n2][0][nCount - 1]) / 2;
                        m_pdZt[i][j][nCount] = 0;
                    }

                    // bottom side
                    else if (i > 0 && i < n1 && j == 0)
                    {
                        if (abs(m_pdYt[i][j][nCount - 1]) <= 0.25 * 0.001)
                        {
                            m_pdXt[i][j][nCount] = m_dAlpha * (m_pdXt[i + 1][j][nCount - 1] + m_pdXt[i][j][nCount - 1] + m_pdXt[i - 1][j][nCount - 1]) / 3 + (1 - m_dAlpha) * (m_pdXt[i][j + 1][nCount - 1] + m_pdZt[i][j + 1][nCount - 1] *
                                cos(Theta(i, j, nCount)) / sin(Theta(i, j, nCount)) * NormalVectorBottomX(i, nCount));
                            m_pdYt[i][j][nCount] = m_dAlpha * (m_pdYt[i + 1][j][nCount - 1] + m_pdYt[i][j][nCount - 1] + m_pdYt[i - 1][j][nCount - 1]) / 3 + (1 - m_dAlpha) * (m_pdYt[i][j + 1][nCount - 1] + m_pdZt[i][j + 1][nCount - 1] *
                                cos(Theta(i, j, nCount)) / sin(Theta(i, j, nCount)) * NormalVectorBottomY(i, nCount));
                            m_pdZt[i][j][nCount] = 0;
                        }
                        else if ((abs(m_pdYt[i][j][nCount - 1]) <= 1.25 * 0.001) && (abs(m_pdYt[i][j][nCount - 1]) >= 0.75 * 0.001))
                        {
                            m_pdXt[i][j][nCount] = m_dAlpha * (m_pdXt[i + 1][j][nCount - 1] + m_pdXt[i][j][nCount - 1] + m_pdXt[i - 1][j][nCount - 1]) / 3 + (1 - m_dAlpha) * (m_pdXt[i][j + 1][nCount - 1] + m_pdZt[i][j + 1][nCount - 1] *
                                cos(Theta(i, j, nCount)) / sin(Theta(i, j, nCount)) * NormalVectorBottomX(i, nCount));
                            m_pdYt[i][j][nCount] = m_dAlpha * (m_pdYt[i + 1][j][nCount - 1] + m_pdYt[i][j][nCount - 1] + m_pdYt[i - 1][j][nCount - 1]) / 3 + (1 - m_dAlpha) * (m_pdYt[i][j + 1][nCount - 1] + m_pdZt[i][j + 1][nCount - 1] *
                                cos(Theta(i, j, nCount)) / sin(Theta(i, j, nCount)) * NormalVectorBottomY(i, nCount));
                            m_pdZt[i][j][nCount] = 0;
                        }
                        else
                        {
                            //new time=old time [pinned]
                            m_pdXt[i][j][nCount] = m_pdXt[i][j][nCount - 1];
                            m_pdYt[i][j][nCount] = m_pdYt[i][j][nCount - 1];
                            m_pdZt[i][j][nCount] = m_pdZt[i][j][nCount - 1];
                        }
                    }


                    // right side
                    else if (i == n1 && j > 0 && j < m1)
                    {
                        if (abs(m_pdYt[i][j][nCount - 1]) <= 0.25 * 0.001)
                        {
                            m_pdXt[i][j][nCount] = m_dAlpha * (m_pdXt[i][j + 1][nCount - 1] + m_pdXt[i][j][nCount - 1] + m_pdXt[i][j - 1][nCount - 1]) / 3 + (1 - m_dAlpha) * (m_pdXt[i - 1][j][nCount - 1] + m_pdZt[i - 1][j][nCount - 1] *
                                cos(Theta(i, j, nCount)) / sin(Theta(i, j, nCount)) * NormalVectorRightX(j, nCount));
                            m_pdYt[i][j][nCount] = m_dAlpha * (m_pdYt[i][j + 1][nCount - 1] + m_pdYt[i][j][nCount - 1] + m_pdYt[i][j - 1][nCount - 1]) / 3 + (1 - m_dAlpha) * (m_pdYt[i - 1][j][nCount - 1] + m_pdZt[i - 1][j][nCount - 1] *
                                cos(Theta(i, j, nCount)) / sin(Theta(i, j, nCount)) * NormalVectorRightY(j, nCount));
                            m_pdZt[i][j][nCount] = 0;
                        }
                        else if ((abs(m_pdYt[i][j][nCount - 1]) <= 1.25 * 0.001) && (abs(m_pdYt[i][j][nCount - 1]) >= 0.75 * 0.001))
                        {
                            m_pdXt[i][j][nCount] = m_dAlpha * (m_pdXt[i][j + 1][nCount - 1] + m_pdXt[i][j][nCount - 1] + m_pdXt[i][j - 1][nCount - 1]) / 3 + (1 - m_dAlpha) * (m_pdXt[i - 1][j][nCount - 1] + m_pdZt[i - 1][j][nCount - 1] *
                                cos(Theta(i, j, nCount)) / sin(Theta(i, j, nCount)) * NormalVectorRightX(j, nCount));
                            m_pdYt[i][j][nCount] = m_dAlpha * (m_pdYt[i][j + 1][nCount - 1] + m_pdYt[i][j][nCount - 1] + m_pdYt[i][j - 1][nCount - 1]) / 3 + (1 - m_dAlpha) * (m_pdYt[i - 1][j][nCount - 1] + m_pdZt[i - 1][j][nCount - 1] *
                                cos(Theta(i, j, nCount)) / sin(Theta(i, j, nCount)) * NormalVectorRightY(j, nCount));
                            m_pdZt[i][j][nCount] = 0;
                        }
                        else
                        {
                            //new time=old time [pinned]
                            m_pdXt[i][j][nCount] = m_pdXt[i][j][nCount - 1];
                            m_pdYt[i][j][nCount] = m_pdYt[i][j][nCount - 1];
                            m_pdZt[i][j][nCount] = m_pdZt[i][j][nCount - 1];
                        }

                    }

                    // top side
                    else if (i > 0 && i < n1 && j == m1)
                    {
                        if (abs(m_pdYt[i][j][nCount - 1]) <= 0.25 * 0.001)
                        {
                            m_pdXt[i][j][nCount] = m_dAlpha * (m_pdXt[i + 1][j][nCount - 1] + m_pdXt[i][j][nCount - 1] + m_pdXt[i - 1][j][nCount - 1]) / 3 + (1 - m_dAlpha) * (m_pdXt[i][j - 1][nCount - 1] + m_pdZt[i][j - 1][nCount - 1] *
                                cos(Theta(i, j, nCount)) / sin(Theta(i, j, nCount)) * NormalVectorTopX(i, nCount));
                            m_pdYt[i][j][nCount] = m_dAlpha * (m_pdYt[i + 1][j][nCount - 1] + m_pdYt[i][j][nCount - 1] + m_pdYt[i - 1][j][nCount - 1]) / 3 + (1 - m_dAlpha) * (m_pdYt[i][j - 1][nCount - 1] + m_pdZt[i][j - 1][nCount - 1] *
                                cos(Theta(i, j, nCount)) / sin(Theta(i, j, nCount)) * NormalVectorTopY(i, nCount));
                            m_pdZt[i][j][nCount] = 0;
                        }
                        else if ((abs(m_pdYt[i][j][nCount - 1]) <= 1.25 * 0.001) && (abs(m_pdYt[i][j][nCount - 1]) >= 0.75 * 0.001))
                        {
                            m_pdXt[i][j][nCount] = m_dAlpha * (m_pdXt[i + 1][j][nCount - 1] + m_pdXt[i][j][nCount - 1] + m_pdXt[i - 1][j][nCount - 1]) / 3 + (1 - m_dAlpha) * (m_pdXt[i][j - 1][nCount - 1] + m_pdZt[i][j - 1][nCount - 1] *
                                cos(Theta(i, j, nCount)) / sin(Theta(i, j, nCount)) * NormalVectorTopX(i, nCount));
                            m_pdYt[i][j][nCount] = m_dAlpha * (m_pdYt[i + 1][j][nCount - 1] + m_pdYt[i][j][nCount - 1] + m_pdYt[i - 1][j][nCount - 1]) / 3 + (1 - m_dAlpha) * (m_pdYt[i][j - 1][nCount - 1] + m_pdZt[i][j - 1][nCount - 1] *
                                cos(Theta(i, j, nCount)) / sin(Theta(i, j, nCount)) * NormalVectorTopY(i, nCount));
                            m_pdZt[i][j][nCount] = 0;
                        }
                        else
                        {
                            //new time=old time [pinned]
                            m_pdXt[i][j][nCount] = m_pdXt[i][j][nCount - 1];
                            m_pdYt[i][j][nCount] = m_pdYt[i][j][nCount - 1];
                            m_pdZt[i][j][nCount] = m_pdZt[i][j][nCount - 1];
                        }

                    }

                    // left side
                    else if (i == 0 && j > 0 && j < m1)
                    {
                        if (abs(m_pdYt[i][j][nCount - 1]) <= 0.25 * 0.001)
                        {
                            m_pdXt[i][j][nCount] = m_dAlpha * (m_pdXt[i][j + 1][nCount - 1] + m_pdXt[i][j][nCount - 1] + m_pdXt[i][j - 1][nCount - 1]) / 3 + (1 - m_dAlpha) * (m_pdXt[i + 1][j][nCount - 1] + m_pdZt[i + 1][j][nCount - 1] *
                                cos(Theta(i, j, nCount)) / sin(Theta(i, j, nCount)) * NormalVectorLeftX(j, nCount));
                            m_pdYt[i][j][nCount] = m_dAlpha * (m_pdYt[i][j + 1][nCount - 1] + m_pdYt[i][j][nCount - 1] + m_pdYt[i][j - 1][nCount - 1]) / 3 + (1 - m_dAlpha) * (m_pdYt[i + 1][j][nCount - 1] + m_pdZt[i + 1][j][nCount - 1] *
                                cos(Theta(i, j, nCount)) / sin(Theta(i, j, nCount)) * NormalVectorLeftY(j, nCount));
                            m_pdZt[i][j][nCount] = 0;
                        }
                        else if ((abs(m_pdYt[i][j][nCount - 1]) <= 1.25 * 0.001) && (abs(m_pdYt[i][j][nCount - 1]) >= 0.75 * 0.001))
                        {
                            m_pdXt[i][j][nCount] = m_dAlpha * (m_pdXt[i][j + 1][nCount - 1] + m_pdXt[i][j][nCount - 1] + m_pdXt[i][j - 1][nCount - 1]) / 3 + (1 - m_dAlpha) * (m_pdXt[i + 1][j][nCount - 1] + m_pdZt[i + 1][j][nCount - 1] *
                                cos(Theta(i, j, nCount)) / sin(Theta(i, j, nCount)) * NormalVectorLeftX(j, nCount));
                            m_pdYt[i][j][nCount] = m_dAlpha * (m_pdYt[i][j + 1][nCount - 1] + m_pdYt[i][j][nCount - 1] + m_pdYt[i][j - 1][nCount - 1]) / 3 + (1 - m_dAlpha) * (m_pdYt[i + 1][j][nCount - 1] + m_pdZt[i + 1][j][nCount - 1] *
                                cos(Theta(i, j, nCount)) / sin(Theta(i, j, nCount)) * NormalVectorLeftY(j, nCount));
                            m_pdZt[i][j][nCount] = 0;
                        }
                        else
                        {
                            //new time=old time [pinned]
                            m_pdXt[i][j][nCount] = m_pdXt[i][j][nCount - 1];
                            m_pdYt[i][j][nCount] = m_pdYt[i][j][nCount - 1];
                            m_pdZt[i][j][nCount] = m_pdZt[i][j][nCount - 1];
                        }

                    }


                    else
                    {
                        m_pdXt[i][j][nCount] = m_dDeltaT / m_dFric *
                            (KnownNormalForceX(i, j, nCount) +
                                KnownTangentialForceX(i, j, nCount)) + m_pdXt[i][j][nCount - 1];
                        m_pdYt[i][j][nCount] = m_dDeltaT / m_dFric *
                            (KnownNormalForceY(i, j, nCount) +
                                KnownTangentialForceY(i, j, nCount)) + m_pdYt[i][j][nCount - 1];
                        m_pdZt[i][j][nCount] = m_dDeltaT / m_dFric *
                            (KnownNormalForceZ(i, j, nCount) +
                                KnownTangentialForceZ(i, j, nCount)) + m_pdZt[i][j][nCount - 1];
                    }
                    // save theta into array

                    if (nCount == m_nTimes - 1)
                    {
                        theta[i][j] = Theta(i, j, nCount);
                    }
                }
            }


        }
        // save pressure into array

        if (nCount == 0)
        {
            pdVol[nCount] = m_dV;
            pdPress[nCount] = m_dP;
        }

        // save volume into array

        else
        {
            pdVol[nCount] = CurrentVolume(nCount);
            //std::cout << "count is " << nCount << " with volume " << pdVol[nCount] << "\n";
            pdPress[nCount] = Pressure(nCount);
        }


    }
}

/******************************************************
    **        The boundary conditions                  **
    ** The normal vectors for all four sides      **
******************************************************/

double Drop::NormalVectorBottomX(int i, int count)
{
    double x;
    // arc length formula
    x = -(m_pdYt[i + 1][0][count - 1] - m_pdYt[i - 1][0][count - 1]) /
        sqrt((m_pdYt[i + 1][0][count - 1] - m_pdYt[i - 1][0][count - 1]) *
            (m_pdYt[i + 1][0][count - 1] - m_pdYt[i - 1][0][count - 1]) +
            (m_pdXt[i + 1][0][count - 1] - m_pdXt[i - 1][0][count - 1]) *
            (m_pdXt[i + 1][0][count - 1] - m_pdXt[i - 1][0][count - 1]));

    return x;
}
double Drop::NormalVectorBottomY(int i, int count)
{
    double y;
    y = (m_pdXt[i + 1][0][count - 1] - m_pdXt[i - 1][0][count - 1]) /
        sqrt((m_pdYt[i + 1][0][count - 1] - m_pdYt[i - 1][0][count - 1]) *
            (m_pdYt[i + 1][0][count - 1] - m_pdYt[i - 1][0][count - 1]) +
            (m_pdXt[i + 1][0][count - 1] - m_pdXt[i - 1][0][count - 1]) *
            (m_pdXt[i + 1][0][count - 1] - m_pdXt[i - 1][0][count - 1]));

    return y;
}
double Drop::NormalVectorBottomZ(int i, int count)
{
    return 0;
}

double Drop::NormalVectorRightX(int j, int count)
{
    double x;
    x = -(m_pdYt[m_nPoints - 1][j + 1][count - 1] - m_pdYt[m_nPoints - 1][j - 1][count - 1]) /
        sqrt((m_pdYt[m_nPoints - 1][j + 1][count - 1] - m_pdYt[m_nPoints - 1][j - 1][count - 1]) *
            (m_pdYt[m_nPoints - 1][j + 1][count - 1] - m_pdYt[m_nPoints - 1][j - 1][count - 1]) +
            (m_pdXt[m_nPoints - 1][j + 1][count - 1] - m_pdXt[m_nPoints - 1][j - 1][count - 1]) *
            (m_pdXt[m_nPoints - 1][j + 1][count - 1] - m_pdXt[m_nPoints - 1][j - 1][count - 1]));

    return x;
}
double Drop::NormalVectorRightY(int j, int count)
{
    double y;
    y = (m_pdXt[m_nPoints - 1][j + 1][count - 1] - m_pdXt[m_nPoints - 1][j - 1][count - 1]) /
        sqrt((m_pdYt[m_nPoints - 1][j + 1][count - 1] - m_pdYt[m_nPoints - 1][j - 1][count - 1]) *
            (m_pdYt[m_nPoints - 1][j + 1][count - 1] - m_pdYt[m_nPoints - 1][j - 1][count - 1]) +
            (m_pdXt[m_nPoints - 1][j + 1][count - 1] - m_pdXt[m_nPoints - 1][j - 1][count - 1]) *
            (m_pdXt[m_nPoints - 1][j + 1][count - 1] - m_pdXt[m_nPoints - 1][j - 1][count - 1]));

    return y;
}
double Drop::NormalVectorRightZ(int j, int count)
{
    return 0;
}

double Drop::NormalVectorTopX(int i, int count)
{
    double x;
    x = (m_pdYt[i + 1][m_mPoints - 1][count - 1] - m_pdYt[i - 1][m_mPoints - 1][count - 1]) /
        sqrt((m_pdYt[i + 1][m_mPoints - 1][count - 1] - m_pdYt[i - 1][m_mPoints - 1][count - 1]) *
            (m_pdYt[i + 1][m_mPoints - 1][count - 1] - m_pdYt[i - 1][m_mPoints - 1][count - 1]) +
            (m_pdXt[i + 1][m_mPoints - 1][count - 1] - m_pdXt[i - 1][m_mPoints - 1][count - 1]) *
            (m_pdXt[i + 1][m_mPoints - 1][count - 1] - m_pdXt[i - 1][m_mPoints - 1][count - 1]));

    return x;
}
double Drop::NormalVectorTopY(int i, int count)
{
    double y;
    y = -(m_pdXt[i + 1][m_mPoints - 1][count - 1] - m_pdXt[i - 1][m_mPoints - 1][count - 1]) /
        sqrt((m_pdYt[i + 1][m_mPoints - 1][count - 1] - m_pdYt[i - 1][m_mPoints - 1][count - 1]) *
            (m_pdYt[i + 1][m_mPoints - 1][count - 1] - m_pdYt[i - 1][m_mPoints - 1][count - 1]) +
            (m_pdXt[i + 1][m_mPoints - 1][count - 1] - m_pdXt[i - 1][m_mPoints - 1][count - 1]) *
            (m_pdXt[i + 1][m_mPoints - 1][count - 1] - m_pdXt[i - 1][m_mPoints - 1][count - 1]));

    return y;
}
double Drop::NormalVectorTopZ(int i, int count)
{
    return 0;
}

double Drop::NormalVectorLeftX(int j, int count)
{
    double x;
    x = (m_pdYt[0][j + 1][count - 1] - m_pdYt[0][j - 1][count - 1]) /
        sqrt((m_pdYt[0][j + 1][count - 1] - m_pdYt[0][j - 1][count - 1]) *
            (m_pdYt[0][j + 1][count - 1] - m_pdYt[0][j - 1][count - 1]) +
            (m_pdXt[0][j + 1][count - 1] - m_pdXt[0][j - 1][count - 1]) *
            (m_pdXt[0][j + 1][count - 1] - m_pdXt[0][j - 1][count - 1]));

    return x;
}
double Drop::NormalVectorLeftY(int j, int count)
{
    double y;
    y = -(m_pdXt[0][j + 1][count - 1] - m_pdXt[0][j - 1][count - 1]) /
        sqrt((m_pdYt[0][j + 1][count - 1] - m_pdYt[0][j - 1][count - 1]) *
            (m_pdYt[0][j + 1][count - 1] - m_pdYt[0][j - 1][count - 1]) +
            (m_pdXt[0][j + 1][count - 1] - m_pdXt[0][j - 1][count - 1]) *
            (m_pdXt[0][j + 1][count - 1] - m_pdXt[0][j - 1][count - 1]));

    return y;
}
double Drop::NormalVectorLeftZ(int j, int count)
{
    return 0;
}

double Drop::Theta(int i, int j, int count)
{
    double theta;
    double point = fabs(m_pdXt[i][j][count - 1]);
    double qoint = fabs(m_pdYt[i][j][count - 1]);
    double b = 1.0e3;
    //roint = sqrt((point) * (point)+qoint * qoint);

    //    Circle of radius 10 * 10^(-3)

    //theta = -PI * 20. / 180 + (110 - 20) * -PI / 180 * ((1 + tanh(b * (roint - 10.e-4))) / 2);
    //     if( roint < 10.e-4 )
    //             theta=-PI*20./180.;
    //     else
    //            theta=-PI*110./180.;


    //    Square with a side of 20*10^(-3)
    /*    if( point < 10.e-4 && qoint < 10.e-4 )
            theta = -PI*20./180.;
         else
            theta = -PI*110./180.;
    */

    //    Hexagon with a smaller radius of 10 * 10^(-3)
    /*    if( qoint < 10.e-4 &&
            qoint < (15 * point / (5*sqrt(3) - 10) )*10.e-4 )
            theta = -PI*20./180.;
         else
            theta = -PI*110./180.;
    */
    theta = -75 * PI / 180;
    return theta;
}

// the X component of all the forces acting on the point
// in the Normal direction

double Drop::KnownNormalForceX(int i, int j, int count)
{
    double f;
    f = PressureForceX(i, j, count) + SigmaForceX(i, j, count);

    return f;
}

// the Y component of all the forces acting on the point
// in the Normal direction

double Drop::KnownNormalForceY(int i, int j, int count)
{
    double f;
    f = PressureForceY(i, j, count) + SigmaForceY(i, j, count);

    return f;
}

// the Z component of all the forces acting on the point
// in the Normal direction

double Drop::KnownNormalForceZ(int i, int j, int count)
{
    double f;
    f = PressureForceZ(i, j, count) + SigmaForceZ(i, j, count);

    return f;
}

// The X component of the tangential force

double Drop::KnownTangentialForceX(int i, int j, int count)
{
    double f;

    f = TangentPartX(i, j, count) -
        (TangentPartX(i, j, count) * NormalVectorX(i, j, count) +
            TangentPartY(i, j, count) * NormalVectorY(i, j, count) +
            TangentPartZ(i, j, count) * NormalVectorZ(i, j, count)) *
        NormalVectorX(i, j, count);

    return f;
}

// The Y component of the tangential force

double Drop::KnownTangentialForceY(int i, int j, int count)
{
    double f;

    f = TangentPartY(i, j, count) -
        (TangentPartX(i, j, count) * NormalVectorX(i, j, count) +
            TangentPartY(i, j, count) * NormalVectorY(i, j, count) +
            TangentPartZ(i, j, count) * NormalVectorZ(i, j, count)) *
        NormalVectorY(i, j, count);

    return f;
}

// The Z component of the tangential force

double Drop::KnownTangentialForceZ(int i, int j, int count)
{
    double f;

    f = TangentPartZ(i, j, count) -
        (TangentPartX(i, j, count) * NormalVectorX(i, j, count) +
            TangentPartY(i, j, count) * NormalVectorY(i, j, count) +
            TangentPartZ(i, j, count) * NormalVectorZ(i, j, count)) *
        NormalVectorZ(i, j, count);

    return f;
}

// calculate the current volume at a particular time step
double Drop::CurrentVolume(int count)
{
    double v_current = 0;
    int c1 = count - 1;
    for (int i = 0; i < n1; i++)
    {
        for (int j = 0; j < m1; j++)
        {
            v_current = (double)(v_current +
                ((m_pdXt[i + 1][j][c1] - m_pdXt[i][j][c1]) * (m_pdYt[i][j + 1][c1] - m_pdYt[i][j][c1]) - (m_pdXt[i][j + 1][c1] - m_pdXt[i][j][c1]) * (m_pdYt[i + 1][j][c1] - m_pdYt[i][j][c1]) +
                    (m_pdXt[i + 1][j + 1][c1] - m_pdXt[i][j + 1][c1]) * (m_pdYt[i + 1][j + 1][c1] - m_pdYt[i + 1][j][c1]) - (m_pdXt[i + 1][j + 1][c1] - m_pdXt[i + 1][j][c1]) * (m_pdYt[i + 1][j + 1][c1] - m_pdYt[i][j + 1][c1]))
                * (m_pdZt[i][j][c1] + m_pdZt[i + 1][j][c1] + m_pdZt[i][j + 1][c1] + m_pdZt[i + 1][j + 1][c1]) / 8);
        }
    }
    //std::cout << "count is " << count << " with volume " << v_current << "\n"; //old one
    return v_current;
}

// adjust pressure due to changes in volume

double Drop::Pressure(int count)
{
    double pr;
    double delta_p;

    pr = m_PressureFactor;

    //    delta_p =  1.e4* ( 1.E-3 - m_pdZt[5][5][count-1]);//(m_dV - CurrentVolume( count ));
    delta_p = m_PressRelax * (m_dV - CurrentVolume(count));
    pr += delta_p;

    m_dP = exp(pr) * m_dSigma * pow(m_dV / (CurrentVolume(count) * CurrentVolume(count)),
        (1. + .1 * ((m_dV / CurrentVolume(count) + CurrentVolume(count) / m_dV) - 2)) / 3.);

    //    m_dP += delta_p;
    m_PressureFactor = pr;
    return m_dP;
}

// The X component of the Pressure Force

double Drop::PressureForceX(int j, int k, int count)
{
    double p;
    p = PressureConst(j, k, count) * NormalVectorX(j, k, count);

    return p;
}

// The Y component of the Pressure Force

double Drop::PressureForceY(int j, int k, int count)
{
    double p;
    p = PressureConst(j, k, count) * NormalVectorY(j, k, count);

    return p;
}

// The Z component of the Pressure Force

double Drop::PressureForceZ(int j, int k, int count)
{
    double p;
    p = PressureConst(j, k, count) * NormalVectorZ(j, k, count);

    return p;
}

// The X component of the Force due to Surface Tention

double Drop::SigmaForceX(int i, int j, int count)
{
    double s;
    s = m_dSigma * (NormalVectorX(i, j, count) * SigmaDiffX(i, j, count) +
        NormalVectorY(i, j, count) * SigmaDiffY(i, j, count) +
        NormalVectorZ(i, j, count) * SigmaDiffZ(i, j, count)) *
        NormalVectorX(i, j, count);


    return s;
}

// The Y component of the Force due to Surface Tention

double Drop::SigmaForceY(int i, int j, int count)
{
    double s;
    s = m_dSigma * (NormalVectorX(i, j, count) * SigmaDiffX(i, j, count) +
        NormalVectorY(i, j, count) * SigmaDiffY(i, j, count) +
        NormalVectorZ(i, j, count) * SigmaDiffZ(i, j, count)) *
        NormalVectorY(i, j, count);


    return s;
}

// The Z component of the Force due to Surface Tention

double Drop::SigmaForceZ(int i, int j, int count)
{
    double s;
    s = m_dSigma * (NormalVectorX(i, j, count) * SigmaDiffX(i, j, count) +
        NormalVectorY(i, j, count) * SigmaDiffY(i, j, count) +
        NormalVectorZ(i, j, count) * SigmaDiffZ(i, j, count)) *
        NormalVectorZ(i, j, count);


    return s;
}

// The X component of the part of the tangential force

double Drop::TangentPartX(int i, int j, int count)
{
    double w;

    w = m_dTan * (
        (m_pdXt[i][j - 1][count - 1] - m_pdXt[i][j][count - 1]) +
        (m_pdXt[i][j + 1][count - 1] - m_pdXt[i][j][count - 1]) +
        (m_pdXt[i + 1][j][count - 1] - m_pdXt[i][j][count - 1]) +
        (m_pdXt[i - 1][j][count - 1] - m_pdXt[i][j][count - 1]));

    return w;
}

// The Y component of the part of the tangential force

double Drop::TangentPartY(int i, int j, int count)
{
    double w;

    w = m_dTan * (
        (m_pdYt[i][j - 1][count - 1] - m_pdYt[i][j][count - 1]) +
        (m_pdYt[i][j + 1][count - 1] - m_pdYt[i][j][count - 1]) +
        (m_pdYt[i + 1][j][count - 1] - m_pdYt[i][j][count - 1]) +
        (m_pdYt[i - 1][j][count - 1] - m_pdYt[i][j][count - 1]));

    return w;
}

// The Z component of the part of the tangential force

double Drop::TangentPartZ(int i, int j, int count)
{
    double w;

    w = m_dTan * (
        (m_pdZt[i][j - 1][count - 1] - m_pdZt[i][j][count - 1]) +
        (m_pdZt[i][j + 1][count - 1] - m_pdZt[i][j][count - 1]) +
        (m_pdZt[i + 1][j][count - 1] - m_pdZt[i][j][count - 1]) +
        (m_pdZt[i - 1][j][count - 1] - m_pdZt[i][j][count - 1]));

    return w;
}

// the X coordinate of the normal vector to the surface at the
// point i,j at time count

double Drop::NormalVectorX(int i, int j, int count)
{
    double n;
    int c1 = count - 1;
    n = ((m_pdYt[i + 1][j][c1] - m_pdYt[i - 1][j][c1]) *
        (m_pdZt[i][j + 1][c1] - m_pdZt[i][j - 1][c1]) -
        (m_pdZt[i + 1][j][c1] - m_pdZt[i - 1][j][c1]) *
        (m_pdYt[i][j + 1][c1] - m_pdYt[i][j - 1][c1]))
        /
        sqrt(((m_pdYt[i + 1][j][c1] - m_pdYt[i - 1][j][c1]) *
            (m_pdZt[i][j + 1][c1] - m_pdZt[i][j - 1][c1]) -
            (m_pdZt[i + 1][j][c1] - m_pdZt[i - 1][j][c1]) *
            (m_pdYt[i][j + 1][c1] - m_pdYt[i][j - 1][c1]))
            *
            ((m_pdYt[i + 1][j][c1] - m_pdYt[i - 1][j][c1]) *
                (m_pdZt[i][j + 1][c1] - m_pdZt[i][j - 1][c1]) -
                (m_pdZt[i + 1][j][c1] - m_pdZt[i - 1][j][c1]) *
                (m_pdYt[i][j + 1][c1] - m_pdYt[i][j - 1][c1]))
            +
            ((m_pdXt[i + 1][j][c1] - m_pdXt[i - 1][j][c1]) *
                (m_pdZt[i][j + 1][c1] - m_pdZt[i][j - 1][c1]) -
                (m_pdZt[i + 1][j][c1] - m_pdZt[i - 1][j][c1]) *
                (m_pdXt[i][j + 1][c1] - m_pdXt[i][j - 1][c1]))
            *
            ((m_pdXt[i + 1][j][c1] - m_pdXt[i - 1][j][c1]) *
                (m_pdZt[i][j + 1][c1] - m_pdZt[i][j - 1][c1]) -
                (m_pdZt[i + 1][j][c1] - m_pdZt[i - 1][j][c1]) *
                (m_pdXt[i][j + 1][c1] - m_pdXt[i][j - 1][c1]))
            +
            ((m_pdXt[i + 1][j][c1] - m_pdXt[i - 1][j][c1]) *
                (m_pdYt[i][j + 1][c1] - m_pdYt[i][j - 1][c1]) -
                (m_pdYt[i + 1][j][c1] - m_pdYt[i - 1][j][c1]) *
                (m_pdXt[i][j + 1][c1] - m_pdXt[i][j - 1][c1]))
            *
            ((m_pdXt[i + 1][j][c1] - m_pdXt[i - 1][j][c1]) *
                (m_pdYt[i][j + 1][c1] - m_pdYt[i][j - 1][c1]) -
                (m_pdYt[i + 1][j][c1] - m_pdYt[i - 1][j][c1]) *
                (m_pdXt[i][j + 1][c1] - m_pdXt[i][j - 1][c1])));

    return n;
}

// the Y coordinate of the normal vector to the surface at the
// point i,j at time count

double Drop::NormalVectorY(int i, int j, int count)
{
    int c1 = count - 1;
    double n;
    n = ((m_pdZt[i + 1][j][c1] - m_pdZt[i - 1][j][c1]) *
        (m_pdXt[i][j + 1][c1] - m_pdXt[i][j - 1][c1]) -
        (m_pdXt[i + 1][j][c1] - m_pdXt[i - 1][j][c1]) *
        (m_pdZt[i][j + 1][c1] - m_pdZt[i][j - 1][c1]))
        /
        sqrt(((m_pdYt[i + 1][j][c1] - m_pdYt[i - 1][j][c1]) *
            (m_pdZt[i][j + 1][c1] - m_pdZt[i][j - 1][c1]) -
            (m_pdZt[i + 1][j][c1] - m_pdZt[i - 1][j][c1]) *
            (m_pdYt[i][j + 1][c1] - m_pdYt[i][j - 1][c1]))
            *
            ((m_pdYt[i + 1][j][c1] - m_pdYt[i - 1][j][c1]) *
                (m_pdZt[i][j + 1][c1] - m_pdZt[i][j - 1][c1]) -
                (m_pdZt[i + 1][j][c1] - m_pdZt[i - 1][j][c1]) *
                (m_pdYt[i][j + 1][c1] - m_pdYt[i][j - 1][c1]))
            +
            ((m_pdXt[i + 1][j][c1] - m_pdXt[i - 1][j][c1]) *
                (m_pdZt[i][j + 1][c1] - m_pdZt[i][j - 1][c1]) -
                (m_pdZt[i + 1][j][c1] - m_pdZt[i - 1][j][c1]) *
                (m_pdXt[i][j + 1][c1] - m_pdXt[i][j - 1][c1]))
            *
            ((m_pdXt[i + 1][j][c1] - m_pdXt[i - 1][j][c1]) *
                (m_pdZt[i][j + 1][c1] - m_pdZt[i][j - 1][c1]) -
                (m_pdZt[i + 1][j][c1] - m_pdZt[i - 1][j][c1]) *
                (m_pdXt[i][j + 1][c1] - m_pdXt[i][j - 1][c1]))
            +
            ((m_pdXt[i + 1][j][c1] - m_pdXt[i - 1][j][c1]) *
                (m_pdYt[i][j + 1][c1] - m_pdYt[i][j - 1][c1]) -
                (m_pdYt[i + 1][j][c1] - m_pdYt[i - 1][j][c1]) *
                (m_pdXt[i][j + 1][c1] - m_pdXt[i][j - 1][c1]))
            *
            ((m_pdXt[i + 1][j][c1] - m_pdXt[i - 1][j][c1]) *
                (m_pdYt[i][j + 1][c1] - m_pdYt[i][j - 1][c1]) -
                (m_pdYt[i + 1][j][c1] - m_pdYt[i - 1][j][c1]) *
                (m_pdXt[i][j + 1][c1] - m_pdXt[i][j - 1][c1])));

    return n;
}

// the Z coordinate of the normal vector to the surface at the
// point i,j at time count

double Drop::NormalVectorZ(int i, int j, int count)
{
    int c1 = count - 1;
    double n;
    n = ((m_pdXt[i + 1][j][c1] - m_pdXt[i - 1][j][c1]) *
        (m_pdYt[i][j + 1][c1] - m_pdYt[i][j - 1][c1]) -
        (m_pdYt[i + 1][j][c1] - m_pdYt[i - 1][j][c1]) *
        (m_pdXt[i][j + 1][c1] - m_pdXt[i][j - 1][c1]))
        /
        sqrt(((m_pdYt[i + 1][j][c1] - m_pdYt[i - 1][j][c1]) *
            (m_pdZt[i][j + 1][c1] - m_pdZt[i][j - 1][c1]) -
            (m_pdZt[i + 1][j][c1] - m_pdZt[i - 1][j][c1]) *
            (m_pdYt[i][j + 1][c1] - m_pdYt[i][j - 1][c1]))
            *
            ((m_pdYt[i + 1][j][c1] - m_pdYt[i - 1][j][c1]) *
                (m_pdZt[i][j + 1][c1] - m_pdZt[i][j - 1][c1]) -
                (m_pdZt[i + 1][j][c1] - m_pdZt[i - 1][j][c1]) *
                (m_pdYt[i][j + 1][c1] - m_pdYt[i][j - 1][c1]))
            +
            ((m_pdXt[i + 1][j][c1] - m_pdXt[i - 1][j][c1]) *
                (m_pdZt[i][j + 1][c1] - m_pdZt[i][j - 1][c1]) -
                (m_pdZt[i + 1][j][c1] - m_pdZt[i - 1][j][c1]) *
                (m_pdXt[i][j + 1][c1] - m_pdXt[i][j - 1][c1]))
            *
            ((m_pdXt[i + 1][j][c1] - m_pdXt[i - 1][j][c1]) *
                (m_pdZt[i][j + 1][c1] - m_pdZt[i][j - 1][c1]) -
                (m_pdZt[i + 1][j][c1] - m_pdZt[i - 1][j][c1]) *
                (m_pdXt[i][j + 1][c1] - m_pdXt[i][j - 1][c1]))
            +
            ((m_pdXt[i + 1][j][c1] - m_pdXt[i - 1][j][c1]) *
                (m_pdYt[i][j + 1][c1] - m_pdYt[i][j - 1][c1]) -
                (m_pdYt[i + 1][j][c1] - m_pdYt[i - 1][j][c1]) *
                (m_pdXt[i][j + 1][c1] - m_pdXt[i][j - 1][c1]))
            *
            ((m_pdXt[i + 1][j][c1] - m_pdXt[i - 1][j][c1]) *
                (m_pdYt[i][j + 1][c1] - m_pdYt[i][j - 1][c1]) -
                (m_pdYt[i + 1][j][c1] - m_pdYt[i - 1][j][c1]) *
                (m_pdXt[i][j + 1][c1] - m_pdXt[i][j - 1][c1])));

    return n;
}

// The scalar component of the Pressure Force
/******************************************************
    **                    FORMULA                         **
    ** P/4 * {|(Vi,j-1 - Vi-1,j)X(Vi-1,j - Vi,j+1 )| **
    ** + |(Vi,j-1 - Vi+1,j)X(Vi+1,j - Vi,j+1 )|}     **
 ******************************************************/


double Drop::PressureConst(int i, int j, int count)
{
    double p;
    p = Pressure(count) / 4 * (
        sqrt(((m_pdYt[i][j - 1][count - 1] - m_pdYt[i - 1][j][count - 1]) *
            (m_pdZt[i - 1][j][count - 1] - m_pdZt[i][j + 1][count - 1]) -
            (m_pdZt[i][j - 1][count - 1] - m_pdZt[i - 1][j][count - 1]) *
            (m_pdYt[i - 1][j][count - 1] - m_pdYt[i][j + 1][count - 1]))
            *
            ((m_pdYt[i][j - 1][count - 1] - m_pdYt[i - 1][j][count - 1]) *
                (m_pdZt[i - 1][j][count - 1] - m_pdZt[i][j + 1][count - 1]) -
                (m_pdZt[i][j - 1][count - 1] - m_pdZt[i - 1][j][count - 1]) *
                (m_pdYt[i - 1][j][count - 1] - m_pdYt[i][j + 1][count - 1]))
            +
            ((m_pdXt[i][j - 1][count - 1] - m_pdXt[i - 1][j][count - 1]) *
                (m_pdZt[i - 1][j][count - 1] - m_pdZt[i][j + 1][count - 1]) -
                (m_pdZt[i][j - 1][count - 1] - m_pdZt[i - 1][j][count - 1]) *
                (m_pdXt[i - 1][j][count - 1] - m_pdXt[i][j + 1][count - 1]))
            *
            ((m_pdXt[i][j - 1][count - 1] - m_pdXt[i - 1][j][count - 1]) *
                (m_pdZt[i - 1][j][count - 1] - m_pdZt[i][j + 1][count - 1]) -
                (m_pdZt[i][j - 1][count - 1] - m_pdZt[i - 1][j][count - 1]) *
                (m_pdXt[i - 1][j][count - 1] - m_pdXt[i][j + 1][count - 1]))
            +
            ((m_pdXt[i][j - 1][count - 1] - m_pdXt[i - 1][j][count - 1]) *
                (m_pdYt[i - 1][j][count - 1] - m_pdYt[i][j + 1][count - 1]) -
                (m_pdYt[i][j - 1][count - 1] - m_pdYt[i - 1][j][count - 1]) *
                (m_pdXt[i - 1][j][count - 1] - m_pdXt[i][j + 1][count - 1]))
            *
            ((m_pdXt[i][j - 1][count - 1] - m_pdXt[i - 1][j][count - 1]) *
                (m_pdYt[i - 1][j][count - 1] - m_pdYt[i][j + 1][count - 1]) -
                (m_pdYt[i][j - 1][count - 1] - m_pdYt[i - 1][j][count - 1]) *
                (m_pdXt[i - 1][j][count - 1] - m_pdXt[i][j + 1][count - 1])))
        +
        sqrt(((m_pdYt[i][j - 1][count - 1] - m_pdYt[i + 1][j][count - 1]) *
            (m_pdZt[i + 1][j][count - 1] - m_pdZt[i][j + 1][count - 1]) -
            (m_pdZt[i][j - 1][count - 1] - m_pdZt[i + 1][j][count - 1]) *
            (m_pdYt[i + 1][j][count - 1] - m_pdYt[i][j + 1][count - 1]))
            *
            ((m_pdYt[i][j - 1][count - 1] - m_pdYt[i + 1][j][count - 1]) *
                (m_pdZt[i + 1][j][count - 1] - m_pdZt[i][j + 1][count - 1]) -
                (m_pdZt[i][j - 1][count - 1] - m_pdZt[i + 1][j][count - 1]) *
                (m_pdYt[i + 1][j][count - 1] - m_pdYt[i][j + 1][count - 1]))
            +
            ((m_pdXt[i][j - 1][count - 1] - m_pdXt[i + 1][j][count - 1]) *
                (m_pdZt[i + 1][j][count - 1] - m_pdZt[i][j + 1][count - 1]) -
                (m_pdZt[i][j - 1][count - 1] - m_pdZt[i + 1][j][count - 1]) *
                (m_pdXt[i + 1][j][count - 1] - m_pdXt[i][j + 1][count - 1]))
            *
            ((m_pdXt[i][j - 1][count - 1] - m_pdXt[i + 1][j][count - 1]) *
                (m_pdZt[i + 1][j][count - 1] - m_pdZt[i][j + 1][count - 1]) -
                (m_pdZt[i][j - 1][count - 1] - m_pdZt[i + 1][j][count - 1]) *
                (m_pdXt[i + 1][j][count - 1] - m_pdXt[i][j + 1][count - 1]))
            +
            ((m_pdYt[i][j - 1][count - 1] - m_pdYt[i + 1][j][count - 1]) *
                (m_pdZt[i + 1][j][count - 1] - m_pdZt[i][j + 1][count - 1]) -
                (m_pdZt[i][j - 1][count - 1] - m_pdZt[i + 1][j][count - 1]) *
                (m_pdYt[i + 1][j][count - 1] - m_pdYt[i][j + 1][count - 1]))
            *
            ((m_pdYt[i][j - 1][count - 1] - m_pdYt[i + 1][j][count - 1]) *
                (m_pdZt[i + 1][j][count - 1] - m_pdZt[i][j + 1][count - 1]) -
                (m_pdZt[i][j - 1][count - 1] - m_pdZt[i + 1][j][count - 1]) *
                (m_pdYt[i + 1][j][count - 1] - m_pdYt[i][j + 1][count - 1]))));

    return p;
}

// The X component of the sum of the four vectors, coming
// from the point ij in the four different directions

double Drop::SigmaDiffX(int i, int j, int count)
{
    double total;
    double left;
    double right;
    double up;
    double down;

    left = sqrt((m_pdXt[i][j + 1][count - 1] - m_pdXt[i][j - 1][count - 1]) *
        (m_pdXt[i][j + 1][count - 1] - m_pdXt[i][j - 1][count - 1])
        +
        (m_pdYt[i][j + 1][count - 1] - m_pdYt[i][j - 1][count - 1]) *
        (m_pdYt[i][j + 1][count - 1] - m_pdYt[i][j - 1][count - 1])
        +
        (m_pdZt[i][j + 1][count - 1] - m_pdZt[i][j - 1][count - 1]) *
        (m_pdZt[i][j + 1][count - 1] - m_pdZt[i][j - 1][count - 1]))
        / 2 *
        (m_pdXt[i - 1][j][count - 1] - m_pdXt[i][j][count - 1]) /
        sqrt((m_pdXt[i - 1][j][count - 1] - m_pdXt[i][j][count - 1]) *
            (m_pdXt[i - 1][j][count - 1] - m_pdXt[i][j][count - 1])
            +
            (m_pdYt[i - 1][j][count - 1] - m_pdYt[i][j][count - 1]) *
            (m_pdYt[i - 1][j][count - 1] - m_pdYt[i][j][count - 1])
            +
            (m_pdZt[i - 1][j][count - 1] - m_pdZt[i][j][count - 1]) *
            (m_pdZt[i - 1][j][count - 1] - m_pdZt[i][j][count - 1]));

    right = sqrt((m_pdXt[i][j + 1][count - 1] - m_pdXt[i][j - 1][count - 1]) *
        (m_pdXt[i][j + 1][count - 1] - m_pdXt[i][j - 1][count - 1])
        +
        (m_pdYt[i][j + 1][count - 1] - m_pdYt[i][j - 1][count - 1]) *
        (m_pdYt[i][j + 1][count - 1] - m_pdYt[i][j - 1][count - 1])
        +
        (m_pdZt[i][j + 1][count - 1] - m_pdZt[i][j - 1][count - 1]) *
        (m_pdZt[i][j + 1][count - 1] - m_pdZt[i][j - 1][count - 1]))
        / 2 *
        (m_pdXt[i + 1][j][count - 1] - m_pdXt[i][j][count - 1]) /
        sqrt((m_pdXt[i + 1][j][count - 1] - m_pdXt[i][j][count - 1]) *
            (m_pdXt[i + 1][j][count - 1] - m_pdXt[i][j][count - 1])
            +
            (m_pdYt[i + 1][j][count - 1] - m_pdYt[i][j][count - 1]) *
            (m_pdYt[i + 1][j][count - 1] - m_pdYt[i][j][count - 1])
            +
            (m_pdZt[i + 1][j][count - 1] - m_pdZt[i][j][count - 1]) *
            (m_pdZt[i + 1][j][count - 1] - m_pdZt[i][j][count - 1]));

    up = sqrt((m_pdXt[i - 1][j][count - 1] - m_pdXt[i + 1][j][count - 1]) *
        (m_pdXt[i - 1][j][count - 1] - m_pdXt[i + 1][j][count - 1])
        +
        (m_pdYt[i - 1][j][count - 1] - m_pdYt[i + 1][j][count - 1]) *
        (m_pdYt[i - 1][j][count - 1] - m_pdYt[i + 1][j][count - 1])
        +
        (m_pdZt[i - 1][j][count - 1] - m_pdZt[i + 1][j][count - 1]) *
        (m_pdZt[i - 1][j][count - 1] - m_pdZt[i + 1][j][count - 1]))
        / 2 *
        (m_pdXt[i][j + 1][count - 1] - m_pdXt[i][j][count - 1]) /
        sqrt((m_pdXt[i][j + 1][count - 1] - m_pdXt[i][j][count - 1]) *
            (m_pdXt[i][j + 1][count - 1] - m_pdXt[i][j][count - 1])
            +
            (m_pdYt[i][j + 1][count - 1] - m_pdYt[i][j][count - 1]) *
            (m_pdYt[i][j + 1][count - 1] - m_pdYt[i][j][count - 1])
            +
            (m_pdZt[i][j + 1][count - 1] - m_pdZt[i][j][count - 1]) *
            (m_pdZt[i][j + 1][count - 1] - m_pdZt[i][j][count - 1]));

    down = sqrt((m_pdXt[i - 1][j][count - 1] - m_pdXt[i + 1][j][count - 1]) *
        (m_pdXt[i - 1][j][count - 1] - m_pdXt[i + 1][j][count - 1])
        +
        (m_pdYt[i - 1][j][count - 1] - m_pdYt[i + 1][j][count - 1]) *
        (m_pdYt[i - 1][j][count - 1] - m_pdYt[i + 1][j][count - 1])
        +
        (m_pdZt[i - 1][j][count - 1] - m_pdZt[i + 1][j][count - 1]) *
        (m_pdZt[i - 1][j][count - 1] - m_pdZt[i + 1][j][count - 1]))
        / 2 *
        (m_pdXt[i][j - 1][count - 1] - m_pdXt[i][j][count - 1]) /
        sqrt((m_pdXt[i][j - 1][count - 1] - m_pdXt[i][j][count - 1]) *
            (m_pdXt[i][j - 1][count - 1] - m_pdXt[i][j][count - 1])
            +
            (m_pdYt[i][j - 1][count - 1] - m_pdYt[i][j][count - 1]) *
            (m_pdYt[i][j - 1][count - 1] - m_pdYt[i][j][count - 1])
            +
            (m_pdZt[i][j - 1][count - 1] - m_pdZt[i][j][count - 1]) *
            (m_pdZt[i][j - 1][count - 1] - m_pdZt[i][j][count - 1]));

    total = left + right + up + down;

    return total;
}

// The Y component of the sum of the four vectors, coming
// from the point ij in the four different directions

double Drop::SigmaDiffY(int i, int j, int count)
{
    double total;
    double left;
    double right;
    double up;
    double down;

    left = sqrt((m_pdXt[i][j + 1][count - 1] - m_pdXt[i][j - 1][count - 1]) *
        (m_pdXt[i][j + 1][count - 1] - m_pdXt[i][j - 1][count - 1])
        +
        (m_pdYt[i][j + 1][count - 1] - m_pdYt[i][j - 1][count - 1]) *
        (m_pdYt[i][j + 1][count - 1] - m_pdYt[i][j - 1][count - 1])
        +
        (m_pdZt[i][j + 1][count - 1] - m_pdZt[i][j - 1][count - 1]) *
        (m_pdZt[i][j + 1][count - 1] - m_pdZt[i][j - 1][count - 1]))
        / 2 *
        (m_pdYt[i - 1][j][count - 1] - m_pdYt[i][j][count - 1]) /
        sqrt((m_pdXt[i - 1][j][count - 1] - m_pdXt[i][j][count - 1]) *
            (m_pdXt[i - 1][j][count - 1] - m_pdXt[i][j][count - 1])
            +
            (m_pdYt[i - 1][j][count - 1] - m_pdYt[i][j][count - 1]) *
            (m_pdYt[i - 1][j][count - 1] - m_pdYt[i][j][count - 1])
            +
            (m_pdZt[i - 1][j][count - 1] - m_pdZt[i][j][count - 1]) *
            (m_pdZt[i - 1][j][count - 1] - m_pdZt[i][j][count - 1]));

    right = sqrt((m_pdXt[i][j + 1][count - 1] - m_pdXt[i][j - 1][count - 1]) *
        (m_pdXt[i][j + 1][count - 1] - m_pdXt[i][j - 1][count - 1])
        +
        (m_pdYt[i][j + 1][count - 1] - m_pdYt[i][j - 1][count - 1]) *
        (m_pdYt[i][j + 1][count - 1] - m_pdYt[i][j - 1][count - 1])
        +
        (m_pdZt[i][j + 1][count - 1] - m_pdZt[i][j - 1][count - 1]) *
        (m_pdZt[i][j + 1][count - 1] - m_pdZt[i][j - 1][count - 1]))
        / 2 *
        (m_pdYt[i + 1][j][count - 1] - m_pdYt[i][j][count - 1]) /
        sqrt((m_pdXt[i + 1][j][count - 1] - m_pdXt[i][j][count - 1]) *
            (m_pdXt[i + 1][j][count - 1] - m_pdXt[i][j][count - 1])
            +
            (m_pdYt[i + 1][j][count - 1] - m_pdYt[i][j][count - 1]) *
            (m_pdYt[i + 1][j][count - 1] - m_pdYt[i][j][count - 1])
            +
            (m_pdZt[i + 1][j][count - 1] - m_pdZt[i][j][count - 1]) *
            (m_pdZt[i + 1][j][count - 1] - m_pdZt[i][j][count - 1]));

    up = sqrt((m_pdXt[i - 1][j][count - 1] - m_pdXt[i + 1][j][count - 1]) *
        (m_pdXt[i - 1][j][count - 1] - m_pdXt[i + 1][j][count - 1])
        +
        (m_pdYt[i - 1][j][count - 1] - m_pdYt[i + 1][j][count - 1]) *
        (m_pdYt[i - 1][j][count - 1] - m_pdYt[i + 1][j][count - 1])
        +
        (m_pdZt[i - 1][j][count - 1] - m_pdZt[i + 1][j][count - 1]) *
        (m_pdZt[i - 1][j][count - 1] - m_pdZt[i + 1][j][count - 1]))
        / 2 *
        (m_pdYt[i][j + 1][count - 1] - m_pdYt[i][j][count - 1]) /
        sqrt((m_pdXt[i][j + 1][count - 1] - m_pdXt[i][j][count - 1]) *
            (m_pdXt[i][j + 1][count - 1] - m_pdXt[i][j][count - 1])
            +
            (m_pdYt[i][j + 1][count - 1] - m_pdYt[i][j][count - 1]) *
            (m_pdYt[i][j + 1][count - 1] - m_pdYt[i][j][count - 1])
            +
            (m_pdZt[i][j + 1][count - 1] - m_pdZt[i][j][count - 1]) *
            (m_pdZt[i][j + 1][count - 1] - m_pdZt[i][j][count - 1]));

    down = sqrt((m_pdXt[i - 1][j][count - 1] - m_pdXt[i + 1][j][count - 1]) *
        (m_pdXt[i - 1][j][count - 1] - m_pdXt[i + 1][j][count - 1])
        +
        (m_pdYt[i - 1][j][count - 1] - m_pdYt[i + 1][j][count - 1]) *
        (m_pdYt[i - 1][j][count - 1] - m_pdYt[i + 1][j][count - 1])
        +
        (m_pdZt[i - 1][j][count - 1] - m_pdZt[i + 1][j][count - 1]) *
        (m_pdZt[i - 1][j][count - 1] - m_pdZt[i + 1][j][count - 1]))
        / 2 *
        (m_pdYt[i][j - 1][count - 1] - m_pdYt[i][j][count - 1]) /
        sqrt((m_pdXt[i][j - 1][count - 1] - m_pdXt[i][j][count - 1]) *
            (m_pdXt[i][j - 1][count - 1] - m_pdXt[i][j][count - 1])
            +
            (m_pdYt[i][j - 1][count - 1] - m_pdYt[i][j][count - 1]) *
            (m_pdYt[i][j - 1][count - 1] - m_pdYt[i][j][count - 1])
            +
            (m_pdZt[i][j - 1][count - 1] - m_pdZt[i][j][count - 1]) *
            (m_pdZt[i][j - 1][count - 1] - m_pdZt[i][j][count - 1]));

    total = left + right + up + down;

    return total;
}

// The Z component of the sum of the four vectors, coming
// from the point ij in the four different directions

double Drop::SigmaDiffZ(int i, int j, int count)
{
    double total;
    double left;
    double right;
    double up;
    double down;

    left = sqrt((m_pdXt[i][j + 1][count - 1] - m_pdXt[i][j - 1][count - 1]) *
        (m_pdXt[i][j + 1][count - 1] - m_pdXt[i][j - 1][count - 1])
        +
        (m_pdYt[i][j + 1][count - 1] - m_pdYt[i][j - 1][count - 1]) *
        (m_pdYt[i][j + 1][count - 1] - m_pdYt[i][j - 1][count - 1])
        +
        (m_pdZt[i][j + 1][count - 1] - m_pdZt[i][j - 1][count - 1]) *
        (m_pdZt[i][j + 1][count - 1] - m_pdZt[i][j - 1][count - 1]))
        / 2 *
        (m_pdZt[i - 1][j][count - 1] - m_pdZt[i][j][count - 1]) /
        sqrt((m_pdXt[i - 1][j][count - 1] - m_pdXt[i][j][count - 1]) *
            (m_pdXt[i - 1][j][count - 1] - m_pdXt[i][j][count - 1])
            +
            (m_pdYt[i - 1][j][count - 1] - m_pdYt[i][j][count - 1]) *
            (m_pdYt[i - 1][j][count - 1] - m_pdYt[i][j][count - 1])
            +
            (m_pdZt[i - 1][j][count - 1] - m_pdZt[i][j][count - 1]) *
            (m_pdZt[i - 1][j][count - 1] - m_pdZt[i][j][count - 1]));

    right = sqrt((m_pdXt[i][j + 1][count - 1] - m_pdXt[i][j - 1][count - 1]) *
        (m_pdXt[i][j + 1][count - 1] - m_pdXt[i][j - 1][count - 1])
        +
        (m_pdYt[i][j + 1][count - 1] - m_pdYt[i][j - 1][count - 1]) *
        (m_pdYt[i][j + 1][count - 1] - m_pdYt[i][j - 1][count - 1])
        +
        (m_pdZt[i][j + 1][count - 1] - m_pdZt[i][j - 1][count - 1]) *
        (m_pdZt[i][j + 1][count - 1] - m_pdZt[i][j - 1][count - 1]))
        / 2 *
        (m_pdZt[i + 1][j][count - 1] - m_pdZt[i][j][count - 1]) /
        sqrt((m_pdXt[i + 1][j][count - 1] - m_pdXt[i][j][count - 1]) *
            (m_pdXt[i + 1][j][count - 1] - m_pdXt[i][j][count - 1])
            +
            (m_pdYt[i + 1][j][count - 1] - m_pdYt[i][j][count - 1]) *
            (m_pdYt[i + 1][j][count - 1] - m_pdYt[i][j][count - 1])
            +
            (m_pdZt[i + 1][j][count - 1] - m_pdZt[i][j][count - 1]) *
            (m_pdZt[i + 1][j][count - 1] - m_pdZt[i][j][count - 1]));

    up = sqrt((m_pdXt[i - 1][j][count - 1] - m_pdXt[i + 1][j][count - 1]) *
        (m_pdXt[i - 1][j][count - 1] - m_pdXt[i + 1][j][count - 1])
        +
        (m_pdYt[i - 1][j][count - 1] - m_pdYt[i + 1][j][count - 1]) *
        (m_pdYt[i - 1][j][count - 1] - m_pdYt[i + 1][j][count - 1])
        +
        (m_pdZt[i - 1][j][count - 1] - m_pdZt[i + 1][j][count - 1]) *
        (m_pdZt[i - 1][j][count - 1] - m_pdZt[i + 1][j][count - 1]))
        / 2 *
        (m_pdZt[i][j + 1][count - 1] - m_pdZt[i][j][count - 1]) /
        sqrt((m_pdXt[i][j + 1][count - 1] - m_pdXt[i][j][count - 1]) *
            (m_pdXt[i][j + 1][count - 1] - m_pdXt[i][j][count - 1])
            +
            (m_pdYt[i][j + 1][count - 1] - m_pdYt[i][j][count - 1]) *
            (m_pdYt[i][j + 1][count - 1] - m_pdYt[i][j][count - 1])
            +
            (m_pdZt[i][j + 1][count - 1] - m_pdZt[i][j][count - 1]) *
            (m_pdZt[i][j + 1][count - 1] - m_pdZt[i][j][count - 1]));

    down = sqrt((m_pdXt[i - 1][j][count - 1] - m_pdXt[i + 1][j][count - 1]) *
        (m_pdXt[i - 1][j][count - 1] - m_pdXt[i + 1][j][count - 1])
        +
        (m_pdYt[i - 1][j][count - 1] - m_pdYt[i + 1][j][count - 1]) *
        (m_pdYt[i - 1][j][count - 1] - m_pdYt[i + 1][j][count - 1])
        +
        (m_pdZt[i - 1][j][count - 1] - m_pdZt[i + 1][j][count - 1]) *
        (m_pdZt[i - 1][j][count - 1] - m_pdZt[i + 1][j][count - 1]))
        / 2 *
        (m_pdZt[i][j - 1][count - 1] - m_pdZt[i][j][count - 1]) /
        sqrt((m_pdXt[i][j - 1][count - 1] - m_pdXt[i][j][count - 1]) *
            (m_pdXt[i][j - 1][count - 1] - m_pdXt[i][j][count - 1])
            +
            (m_pdYt[i][j - 1][count - 1] - m_pdYt[i][j][count - 1]) *
            (m_pdYt[i][j - 1][count - 1] - m_pdYt[i][j][count - 1])
            +
            (m_pdZt[i][j - 1][count - 1] - m_pdZt[i][j][count - 1]) *
            (m_pdZt[i][j - 1][count - 1] - m_pdZt[i][j][count - 1]));

    total = left + right + up + down;

    return total;
}

void Drop::Initialization()
{
    // allocate the memory for the time arrays
    // we can potentially have 2 arrays of arrays of size MAX_POINTS.
    // however, we specify a number of points we are going to have,
    // which is m_nPoints.
    // each element in 2 arrays is an array of size m_nTimes.
    // therefore, we end up having 2*m_nPoints arrays of size m_nTimes.
    // we will store the values of one point at all the times delta t.

    for (int i = 0; i < m_nPoints; i++)
    {
        for (int j = 0; j < m_mPoints; j++)
        {
            m_pdXt[i][j] = new double[m_nTimes];
            m_pdYt[i][j] = new double[m_nTimes];
            m_pdZt[i][j] = new double[m_nTimes];
            for (int k = 0; k < m_nTimes; k++)
            {
                m_pdXt[i][j][k] = 0;
                m_pdYt[i][j][k] = 0;
                m_pdZt[i][j][k] = 0;
            }
        }
    }


    // make arrays initialization
    // set the initial coordinates of the masses.
    // positions at delta time = 0.
    for (int k = 0; k < m_nPoints; k++)
    {
        for (int j = 0; j < m_mPoints; j++)
        {
            int nCount = 0;
            if (k == 0)
                m_pdXt[k][j][nCount] = -0.5;
            else
                m_pdXt[k][j][nCount] = m_pdXt[k - 1][j][nCount] + 1.0 / (m_nPoints - 1);

            if (j == 0)
                m_pdYt[k][j][nCount] = -0.5;
            else
                m_pdYt[k][j][nCount] = m_pdYt[k][j - 1][nCount] + 1.0 / (m_mPoints - 1);

            m_pdZt[k][j][nCount] = m_dV * 36 *
                (m_pdXt[k][j][nCount] * m_pdXt[k][j][nCount] - 0.25) *
                (m_pdYt[k][j][nCount] * m_pdYt[k][j][nCount] - 0.25);
        }
    }

}

void Drop::printCurrentMassInformation(std::string s)
{
    std::string end = ".txt";
    std::ofstream myfile;
    std::string beg = "";
    beg = s;

    for (int count = 0; count < m_nTimes; count = count + 100) //for (int count = 0; count < m_nTimes; count = count + 100)
    {
        beg = s;
        myfile.open(beg.append("X").append(std::to_string(count)).append(end));

        for (int i = 0; i < m_nPoints; i++)
        {
            for (int j = 0; j < m_mPoints; j++)
            {
                if (j == m_mPoints - 1)
                    myfile << m_pdXt[i][j][count] << "\n";
                else
                    myfile << m_pdXt[i][j][count] << " ";
            }
        }

        myfile.close();
        beg = s;
        myfile.open(beg.append("Y").append(std::to_string(count)).append(end));


        for (int i = 0; i < m_nPoints; i++)
        {
            for (int j = 0; j < m_mPoints; j++)
            {
                if (j == m_mPoints - 1)
                    myfile << m_pdYt[i][j][count] << "\n";
                else
                    myfile << m_pdYt[i][j][count] << " ";
            }
        }

        myfile.close();
        beg = s;
        myfile.open(beg.append("Z").append(std::to_string(count)).append(end));

        for (int i = 0; i < m_nPoints; i++)
        {
            for (int j = 0; j < m_mPoints; j++)
            {
                if (j == m_mPoints - 1)
                    myfile << m_pdZt[i][j][count] << "\n";
                else
                    myfile << m_pdZt[i][j][count] << " ";
            }
        }

        myfile.close();
    }

}

Drop::Drop() {

    // TODO: add construction code here


    // number of points (masses) to put on the curve
    m_nPoints = 41; // odd number greater than 7
    m_mPoints = 41; // odd number greater than 7
    n1 = m_nPoints - 1;
    n2 = m_nPoints - 2;
    n3 = m_nPoints - 3;
    m1 = m_mPoints - 1;
    m2 = m_mPoints - 2;
    m3 = m_mPoints - 3;

    // the arrays of X and Y components of coordinates of masses
    // on the curve at the given time.
    // these are the arrays that are being graphed.


    m_dXmin = -5 * 1.E-3;//1.E-1;
    m_dXmax = 5 * 1.E-3;//1.E-1;
    m_dYmin = -1 * 1.E-5;
    m_dYmax = 5 * 1.E-3;

    m_dV = 3E-9;//6.79347E-8;//1.43014E-8; //5.0E-9;        // specify a volume of a drop in picoliters (10^-12)
    m_dTheta = -PI * 75 / 180;        // specify an angle

        // number of delta t's we have
    m_nTimes = 2101;//1101;//1100; //1101 originally
    m_dDeltaT = 5.E-8;
    m_dTan = 1 * 1.E2;
    m_dFric = 4.e-5;
    m_dMass = 0;
    m_dSigma = 60;
    m_dSigma = 60;
    m_dP = 1.E5;
    m_PressureFactor = 0;
    m_PressRelax = 1.0e4;
    m_dAlpha = 0.5;

    for (int i = 0; i < MAX_POINTS; i++)
    {
        for (int j = 0; j < MAX_POINTS; j++)
        {
            m_pdXt[i][j] = 0;
            m_pdYt[i][j] = 0;
            m_pdZt[i][j] = 0;

            theta[i][j] = 0;
        }
    }

    for (int j = 0; j < MAX_TIMES; j++)
    {
        pdPress[j] = 0;
        pdVol[j] = 0;
    }
}
