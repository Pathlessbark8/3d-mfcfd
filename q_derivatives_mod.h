#pragma once
#include "data_structure_mod.h"
#include "q_variables_mod.h"

//
//
void eval_q_derivatives()
//
//
{
    //
    //
    int i, j, k, r, nbh;
    double x_i, y_i, z_i, x_k, y_k, z_k;
    double delx, dely, delz, dist, weights;
    double sum_delx_sqr, sum_dely_sqr, sum_delz_sqr;
    double sum_delx_dely, sum_dely_delz, sum_delz_delx;
    double sum_delx_delq[5] = {0}, sum_dely_delq[5] = {0}, sum_delz_delq[5] = {0};
    double det, one_by_det;
    double temp[5], qtilde_i[5], qtilde_nbh[5];
    //
    for (i = 0; i < max_points; i++)
    //
    {
        x_i = point.x[i];
        y_i = point.y[i];
        z_i = point.z[i];
        //
        sum_delx_sqr = 0.0;
        sum_dely_sqr = 0.0;
        sum_delz_sqr = 0.0;
        //
        sum_delx_dely = 0.0;
        sum_dely_delz = 0.0;
        sum_delz_delx = 0.0;

        for (r = 0; r < 5; r++)
        {
            sum_delx_delq[r] = 0;
            sum_dely_delq[r] = 0;
            sum_delz_delq[r] = 0;
        }
        //        for (r = 0; r < 5; r++)
        //         {
        //             cout << sum_delx_delq[r] << " ";
        //         }
        //         cout << endl;
        //Commented out as array needs to be initialised as 0 which is done when the array is declared
        // sum_delx_delq = 0.0;
        // sum_dely_delq = 0.0;
        // sum_delz_delq = 0.0;
        //
        for (k = 0; k < 5; k++)
        {
            point.qm[0][k][i] = point.q[k][i]; // q_maximum ..
            point.qm[1][k][i] = point.q[k][i]; // q_minimum ..
        }
        for (k = 0; k < point.nbhs[i]; k++)
        //
        {
            nbh = point.conn[k][i];
            //
            for (r = 0; r < 5; r++)
            {
                if (point.q[r][nbh] > point.qm[0][r][i])
                {
                    point.qm[0][r][i] = point.q[r][nbh];
                }
                if (point.q[r][nbh] < point.qm[1][r][i])
                {
                    point.qm[1][r][i] = point.q[r][nbh];
                }
            }
            //
            x_k = point.x[nbh];
            y_k = point.y[nbh];
            z_k = point.z[nbh];
            //
            delx = x_k - x_i;
            dely = y_k - y_i;
            delz = z_k - z_i;
            //
            dist = sqrt(delx * delx + dely * dely + delz * delz);
            weights = 1.00 / (pow(dist, power));
            //          if(i==4)
            //         {cout<<delx<<" "<<dely<<" "<<delz;
            // cout<<endl;}
            //
            sum_delx_sqr = sum_delx_sqr + delx * delx * weights;
            sum_dely_sqr = sum_dely_sqr + dely * dely * weights;
            sum_delz_sqr = sum_delz_sqr + delz * delz * weights;
            //
            sum_delx_dely = sum_delx_dely + delx * dely * weights;
            sum_dely_delz = sum_dely_delz + dely * delz * weights;
            sum_delz_delx = sum_delz_delx + delz * delx * weights;
            //
            for (r = 0; r < 5; r++)
            {
                temp[r] = (point.q[r][nbh] - point.q[r][i]);

                sum_delx_delq[r] = sum_delx_delq[r] + weights * delx * temp[r];
                sum_dely_delq[r] = sum_dely_delq[r] + weights * dely * temp[r];
                sum_delz_delq[r] = sum_delz_delq[r] + weights * delz * temp[r];
            }
            // if (k == 4)
            // {
            //     for (r = 0; r < 5; r++)
            //     {
            //         cout << sum_delx_delq[r] << " ";
            //     }
            //     cout << endl;
            // }
            //
        }
        //
        det = sum_delx_sqr * (sum_dely_sqr * sum_delz_sqr - sum_dely_delz * sum_dely_delz) - sum_delx_dely * (sum_delx_dely * sum_delz_sqr - sum_dely_delz * sum_delz_delx) + sum_delz_delx * (sum_delx_dely * sum_dely_delz - sum_dely_sqr * sum_delz_delx);
        //
        one_by_det = 1.00 / det;
        //
        // cout<<det<<endl;
        //
        for (k = 0; k < 5; k++)
        {
            temp[k] = sum_delx_delq[k] * (sum_dely_sqr * sum_delz_sqr - sum_dely_delz * sum_dely_delz) - sum_dely_delq[k] * (sum_delx_dely * sum_delz_sqr - sum_delz_delx * sum_dely_delz) + sum_delz_delq[k] * (sum_delx_dely * sum_dely_delz - sum_delz_delx * sum_dely_sqr);
        }
        // for (k = 0; k < 5; k++)
        // {
        //     cout<<sum_delx_delq[k] <<" ";
        // }
        // cout<<endl;
        //
        for (k = 0; k < 5; k++)
        {
            point.dq[0][k][i] = temp[k] * one_by_det;
        }
        //
        for (k = 0; k < 5; k++)
        {
            temp[k] = sum_delx_sqr * (sum_dely_delq[k] * sum_delz_sqr - sum_dely_delz * sum_delz_delq[k]) - sum_delx_dely * (sum_delx_delq[k] * sum_delz_sqr - sum_delz_delx * sum_delz_delq[k]) + sum_delz_delx * (sum_delx_delq[k] * sum_dely_delz - sum_delz_delx * sum_dely_delq[k]);
        }
        //
            //  for (r = 0; r < 5; r++)
            //     {
            //         cout << temp[r] << " ";
            //     }
            //     cout << endl;
        for (k = 0; k < 5; k++)
        {
            point.dq[1][k][i] = temp[k] * one_by_det;
        }
        //
        for (k = 0; k < 5; k++)
        {
            temp[k] = sum_delx_sqr * (sum_dely_sqr * sum_delz_delq[k] - sum_dely_delq[k] * sum_dely_delz) - sum_delx_dely * (sum_delx_dely * sum_delz_delq[k] - sum_delx_delq[k] * sum_dely_delz) + sum_delz_delx * (sum_delx_dely * sum_dely_delq[k] - sum_delx_delq[k] * sum_dely_sqr);
        }
        for (k = 0; k < 5; k++)
        {
            point.dq[2][k][i] = temp[k] * one_by_det;
        }
        //
        //
    }
    //
    //
    //	Inner iterations to compute second order accurate q-derivatives ..
    //
    //
    for (j = 0; j < inner_iterations; j++)
    //
    {
        for (i = 0; i < max_points; i++)
        //
        {
            x_i = point.x[i];
            y_i = point.y[i];
            z_i = point.z[i];
            //
            sum_delx_sqr = 0.0;
            sum_dely_sqr = 0.0;
            sum_delz_sqr = 0.0;
            //
            sum_delx_dely = 0.0;
            sum_dely_delz = 0.0;
            sum_delz_delx = 0.0;
            //
            for (k = 0; k < 5; k++)
            {
                sum_delx_delq[k] = 0.0;
                sum_dely_delq[k] = 0.0;
                sum_delz_delq[k] = 0.0;
            }
            //
            for (k = 0; k < point.nbhs[i]; k++)
            //
            {
                nbh = point.conn[k][i];
                //
                x_k = point.x[nbh];
                y_k = point.y[nbh];
                z_k = point.z[nbh];
                //
                delx = x_k - x_i;
                dely = y_k - y_i;
                delz = z_k - z_i;
                //
                dist = sqrt(delx * delx + dely * dely + delz * delz);
                weights = 1.00 / (pow(dist, power));
                //
                sum_delx_sqr = sum_delx_sqr + delx * delx * weights;
                sum_dely_sqr = sum_dely_sqr + dely * dely * weights;
                sum_delz_sqr = sum_delz_sqr + delz * delz * weights;
                //
                sum_delx_dely = sum_delx_dely + delx * dely * weights;
                sum_dely_delz = sum_dely_delz + dely * delz * weights;
                sum_delz_delx = sum_delz_delx + delz * delx * weights;
                //
                for (r = 0; r < 5; r++)
                {
                    qtilde_i[r] = point.q[r][i] - 0.50 * (delx * point.dq[0][r][i] + dely * point.dq[1][r][i] + delz * point.dq[2][r][i]);
                    qtilde_nbh[r] = point.q[r][nbh] - 0.50 * (delx * point.dq[0][r][nbh] + dely * point.dq[1][r][nbh] + delz * point.dq[2][r][nbh]);
                    temp[r] = qtilde_nbh[r] - qtilde_i[r];
                }
                //
                for (r = 0; r < 5; r++)
                {
                    sum_delx_delq[r] = sum_delx_delq[r] + weights * delx * temp[r];
                    sum_dely_delq[r] = sum_dely_delq[r] + weights * dely * temp[r];
                    sum_delz_delq[r] = sum_delz_delq[r] + weights * delz * temp[r];
                }
                //
            }
            //

            det = sum_delx_sqr * (sum_dely_sqr * sum_delz_sqr - sum_dely_delz * sum_dely_delz) - sum_delx_dely * (sum_delx_dely * sum_delz_sqr - sum_dely_delz * sum_delz_delx) + sum_delz_delx * (sum_delx_dely * sum_dely_delz - sum_dely_sqr * sum_delz_delx);
            //

            one_by_det = 1.00 / det;
            //
            for (r = 0; r < 5; r++)
            {
                temp[r] = sum_delx_delq[r] * (sum_dely_sqr * sum_delz_sqr - sum_dely_delz * sum_dely_delz) - sum_dely_delq[r] * (sum_delx_dely * sum_delz_sqr - sum_delz_delx * sum_dely_delz) + sum_delz_delq[r] * (sum_delx_dely * sum_dely_delz - sum_delz_delx * sum_dely_sqr);
            }
            //
            for (r = 0; r < 5; r++)
            {
                point.temp[0][r][i] = temp[r] * one_by_det;
                //
            }
            for (r = 0; r < 5; r++)
            {
                temp[r] = sum_delx_sqr * (sum_dely_delq[r] * sum_delz_sqr - sum_dely_delz * sum_delz_delq[r]) - sum_delx_dely * (sum_delx_delq[r] * sum_delz_sqr - sum_delz_delx * sum_delz_delq[r]) + sum_delz_delx * (sum_delx_delq[r] * sum_dely_delz - sum_delz_delx * sum_dely_delq[r]);
            }
            //
            for (r = 0; r < 5; r++)
            {
                point.temp[1][r][i] = temp[r] * one_by_det;
            }
            //
            for (r = 0; r < 5; r++)
            {
                temp[r] = sum_delx_sqr * (sum_dely_sqr * sum_delz_delq[r] - sum_dely_delq[r] * sum_dely_delz) - sum_delx_dely * (sum_delx_dely * sum_delz_delq[r] - sum_delx_delq[r] * sum_dely_delz) + sum_delz_delx * (sum_delx_dely * sum_dely_delq[r] - sum_delx_delq[r] * sum_dely_sqr);
            }
            //
            for (r = 0; r < 5; r++)
            {
                point.temp[2][r][i] = temp[r] * one_by_det;
            }
            //
        } //	i loop (for points)..
          //
        for (i = 0; i < max_points; i++)
        {
            for (r = 0; r < 5; r++)
            {
                point.dq[0][r][i] = point.temp[0][r][i];
                point.dq[1][r][i] = point.temp[1][r][i];
                point.dq[2][r][i] = point.temp[2][r][i];
            }
        }
        //
    } //	j loop (for inner iterations)..
      //
      //
}
