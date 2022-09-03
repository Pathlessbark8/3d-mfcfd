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
    int i, k, r, nbh;
    double x_i, y_i, z_i, x_k, y_k, z_k;
    double delx, dely, delz, dist, weights;
    double sum_delx_sqr, sum_dely_sqr, sum_delz_sqr;
    double sum_delx_dely, sum_dely_delz, sum_delz_delx;
    double sum_delx_delq[5] = {0}, sum_dely_delq[5] = {0}, sum_delz_delq[5] = {0};
    double det, one_by_det;
    double temp[5];
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
            //
        }
        //
        det = sum_delx_sqr * (sum_dely_sqr * sum_delz_sqr - sum_dely_delz * sum_dely_delz) - sum_delx_dely * (sum_delx_dely * sum_delz_sqr - sum_dely_delz * sum_delz_delx) + sum_delz_delx * (sum_delx_dely * sum_dely_delz - sum_dely_sqr * sum_delz_delx);
        //
        one_by_det = 1.00 / det;
        //
        for (k = 0; k < 5; k++)
        {
            temp[k] = sum_delx_delq[k] * (sum_dely_sqr * sum_delz_sqr - sum_dely_delz * sum_dely_delz) - sum_dely_delq[k] * (sum_delx_dely * sum_delz_sqr - sum_delz_delx * sum_dely_delz) + sum_delz_delq[k] * (sum_delx_dely * sum_dely_delz - sum_delz_delx * sum_dely_sqr);
        }
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
}

void q_inner_loop()
{
    int i, k, r, nbh;
    double x_i, y_i, z_i, x_k, y_k, z_k;
    double delx, dely, delz, dist, weights;
    double sum_delx_sqr, sum_dely_sqr, sum_delz_sqr;
    double sum_delx_dely, sum_dely_delz, sum_delz_delx;
    double sum_delx_delq[5] = {0}, sum_dely_delq[5] = {0}, sum_delz_delq[5] = {0};
    double det, one_by_det;
    double temp[5], qtilde_i[5], qtilde_nbh[5];
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
    } //
}
void update_inner_loop()
{
    int i, r, j;
    for (j = 0; j < inner_iterations; j++)
    //
    {
        for (i = 0; i < max_points; i++)
        {
            for (r = 0; r < 5; r++)
            {
                point.dq[0][r][i] = point.temp[0][r][i];
                point.dq[1][r][i] = point.temp[1][r][i];
                point.dq[2][r][i] = point.temp[2][r][i];
            }
        }
    }
}

// Cuda function

__global__ void eval_q_derivatives_cuda(points &point, double power)
//
//
{
    //
    //
    int  k, r, nbh;
    double x_i, y_i, z_i, x_k, y_k, z_k;
    double delx, dely, delz, dist, weights;
    double sum_delx_sqr, sum_dely_sqr, sum_delz_sqr;
    double sum_delx_dely, sum_dely_delz, sum_delz_delx;
    double sum_delx_delq[5] = {0}, sum_dely_delq[5] = {0}, sum_delz_delq[5] = {0};
    double det, one_by_det;
    double temp[5];
    //
    int bx = blockIdx.x;
    int tx = threadIdx.x;
    int i = bx * blockDim.x + tx;
    if (i < 0 || i >= max_points)
    {
        return;
    }
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
        //
    }
    //
    det = sum_delx_sqr * (sum_dely_sqr * sum_delz_sqr - sum_dely_delz * sum_dely_delz) - sum_delx_dely * (sum_delx_dely * sum_delz_sqr - sum_dely_delz * sum_delz_delx) + sum_delz_delx * (sum_delx_dely * sum_dely_delz - sum_dely_sqr * sum_delz_delx);
    //
    one_by_det = 1.00 / det;
    //
    for (k = 0; k < 5; k++)
    {
        temp[k] = sum_delx_delq[k] * (sum_dely_sqr * sum_delz_sqr - sum_dely_delz * sum_dely_delz) - sum_dely_delq[k] * (sum_delx_dely * sum_delz_sqr - sum_delz_delx * sum_dely_delz) + sum_delz_delq[k] * (sum_delx_dely * sum_dely_delz - sum_delz_delx * sum_dely_sqr);
    }
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

__global__ void q_inner_loop_cuda(points &point,int power)
{
    int  k, r, nbh;
    double x_i, y_i, z_i, x_k, y_k, z_k;
    double delx, dely, delz, dist, weights;
    double sum_delx_sqr, sum_dely_sqr, sum_delz_sqr;
    double sum_delx_dely, sum_dely_delz, sum_delz_delx;
    double sum_delx_delq[5] = {0}, sum_dely_delq[5] = {0}, sum_delz_delq[5] = {0};
    double det, one_by_det;
    double temp[5], qtilde_i[5], qtilde_nbh[5];

    int bx = blockIdx.x;
    int tx = threadIdx.x;
    int i = bx * blockDim.x + tx;
    if (i < 0 || i >= max_points)
    {
        return;
    }

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
}

__global__ void update_inner_loop_cuda(points &point)
{
    int bx = blockIdx.x;
    int tx = threadIdx.x;
    int i = bx * blockDim.x + tx;
    int r;
    if (i < 0 || i >= max_points)
    {
        return;
    }
    for (r = 0; r < 5; r++)
    {
        point.dq[0][r][i] = point.temp[0][r][i];
        point.dq[1][r][i] = point.temp[1][r][i];
        point.dq[2][r][i] = point.temp[2][r][i];
    }
}

//Multi-node version of Q_DERIVATIVES function

__global__ void eval_q_derivatives_multi_nccl(int myRank,splitPoints *splitPoint, double power,int max_points_on_device,int *globalToLocalIndex,int **globalToGhostIndex,transferPoints **receiveBuffer,int *partVector,transferPoints **sendBuffer)
//
//
{
    //
    //
    int  k, r, nbh;
    double x_i, y_i, z_i, x_k, y_k, z_k;
    double delx, dely, delz, dist, weights;
    double sum_delx_sqr, sum_dely_sqr, sum_delz_sqr;
    double sum_delx_dely, sum_dely_delz, sum_delz_delx;
    double sum_delx_delq[5] = {0}, sum_dely_delq[5] = {0}, sum_delz_delq[5] = {0};
    double det, one_by_det;
    double temp[5];
    //
    int bx = blockIdx.x;
    int tx = threadIdx.x;
    int i = bx * blockDim.x + tx;
    if (i < 0 || i >= max_points_on_device)
    {
        return;
    }
    x_i = splitPoint[i].x;
    y_i = splitPoint[i].y;
    z_i = splitPoint[i].z;
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
    //Commented out as array needs to be initialised as 0 which is done when the array is declared
    // sum_delx_delq = 0.0;
    // sum_dely_delq = 0.0;
    // sum_delz_delq = 0.0;
    //
    for (k = 0; k < 5; k++)
    {
        splitPoint[i].qm[0][k] = splitPoint[i].q[k]; // q_maximum ..
        splitPoint[i].qm[1][k] = splitPoint[i].q[k]; // q_minimum ..
    }

        for (k = 0; k < splitPoint[i].numberOfLocalNbhs ; k++)
        //
        {
            int globalIndex = splitPoint[i].localNbhs[k];
            nbh=globalToLocalIndex[globalIndex];
            
            for (r = 0; r < 5; r++)
            {
                if (splitPoint[nbh].q[r] > splitPoint[i].qm[0][r])
                {
                    splitPoint[i].qm[0][r] = splitPoint[nbh].q[r];
                }
                if (splitPoint[nbh].q[r] < splitPoint[i].qm[1][r])
                {
                    splitPoint[i].qm[1][r] = splitPoint[nbh].q[r];
                }
            }
            //
            x_k = splitPoint[nbh].x;
            y_k = splitPoint[nbh].y;
            z_k = splitPoint[nbh].z;
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
                temp[r] = (splitPoint[nbh].q[r] - splitPoint[i].q[r]);

                sum_delx_delq[r] = sum_delx_delq[r] + weights * delx * temp[r];
                sum_dely_delq[r] = sum_dely_delq[r] + weights * dely * temp[r];
                sum_delz_delq[r] = sum_delz_delq[r] + weights * delz * temp[r];
            }
            
        }

        for (k = 0; k < splitPoint[i].numberOfGhostNbhs ; k++)
        
        {
            nbh = splitPoint[i].ghostNbhs[k];
            int device=partVector[nbh];
            int ghostIndex=globalToGhostIndex[device][nbh];
            //
            for (r = 0; r < 5; r++)
            {
                if (receiveBuffer[device][ghostIndex].q[r] > splitPoint[i].qm[0][r])
                {
                    splitPoint[i].qm[0][r] = receiveBuffer[device][ghostIndex].q[r];
                }
                if (receiveBuffer[device][ghostIndex].q[r] < splitPoint[i].qm[1][r])
                {
                    splitPoint[i].qm[1][r] = receiveBuffer[device][ghostIndex].q[r];
                }
            }
            
            x_k = receiveBuffer[device][ghostIndex].x;
            y_k = receiveBuffer[device][ghostIndex].y;
            z_k = receiveBuffer[device][ghostIndex].z;
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
                temp[r] = (receiveBuffer[device][ghostIndex].q[r] - splitPoint[i].q[r]);
                sum_delx_delq[r] = sum_delx_delq[r] + weights * delx * temp[r];
                sum_dely_delq[r] = sum_dely_delq[r] + weights * dely * temp[r];
                sum_delz_delq[r] = sum_delz_delq[r] + weights * delz * temp[r];
            }
            
        }
        //
        det = sum_delx_sqr * (sum_dely_sqr * sum_delz_sqr - sum_dely_delz * sum_dely_delz) - sum_delx_dely * (sum_delx_dely * sum_delz_sqr - sum_dely_delz * sum_delz_delx) + sum_delz_delx * (sum_delx_dely * sum_dely_delz - sum_dely_sqr * sum_delz_delx);
        //
        one_by_det = 1.00 / det;
        //
        for (k = 0; k < 5; k++)
        {
            temp[k]= sum_delx_delq[k] * (sum_dely_sqr * sum_delz_sqr - sum_dely_delz * sum_dely_delz) - sum_dely_delq[k] * (sum_delx_dely * sum_delz_sqr - sum_delz_delx * sum_dely_delz) + sum_delz_delq[k] * (sum_delx_dely * sum_dely_delz - sum_delz_delx * sum_dely_sqr);
        }
        
        for (k = 0; k < 5; k++)
        {
            splitPoint[i].dq[0][k]=temp[k] *one_by_det;
        }
        
        for (k = 0; k < 5; k++)
        {
            temp[k] = sum_delx_sqr * (sum_dely_delq[k] * sum_delz_sqr - sum_dely_delz * sum_delz_delq[k]) - sum_delx_dely * (sum_delx_delq[k] * sum_delz_sqr - sum_delz_delx * sum_delz_delq[k]) + sum_delz_delx * (sum_delx_delq[k] * sum_dely_delz - sum_delz_delx * sum_dely_delq[k]);
        }
        //
        for (k = 0; k < 5; k++)
        {
            splitPoint[i].dq[1][k] = temp[k] * one_by_det;
        }
        //
        for (k = 0; k < 5; k++)
        {
            temp[k] = sum_delx_sqr * (sum_dely_sqr * sum_delz_delq[k] - sum_dely_delq[k] * sum_dely_delz) - sum_delx_dely * (sum_delx_dely * sum_delz_delq[k] - sum_delx_delq[k] * sum_dely_delz) + sum_delz_delx * (sum_delx_dely * sum_dely_delq[k] - sum_delx_delq[k] * sum_dely_sqr);
        }
        for (k = 0; k < 5; k++)
        {
            splitPoint[i].dq[2][k] = temp[k] * one_by_det;
        }
        //
        if(splitPoint[i].isGhost){
			for(int t=0;t<splitPoint[i].numberOfPartitionsToSendTo;++t){
				for(int r=0;r<5;++r){
					sendBuffer[splitPoint[i].partitions[t]][splitPoint[i].ghostIndex[t]].dq[0][r]=splitPoint[i].dq[0][r];
                    sendBuffer[splitPoint[i].partitions[t]][splitPoint[i].ghostIndex[t]].dq[1][r]=splitPoint[i].dq[1][r];
					sendBuffer[splitPoint[i].partitions[t]][splitPoint[i].ghostIndex[t]].dq[2][r]=splitPoint[i].dq[2][r];
                    sendBuffer[splitPoint[i].partitions[t]][splitPoint[i].ghostIndex[t]].qm[0][r]=splitPoint[i].qm[0][r];
                    sendBuffer[splitPoint[i].partitions[t]][splitPoint[i].ghostIndex[t]].qm[1][r]=splitPoint[i].qm[1][r];
				}
			}
		}
}

__global__ void q_inner_loop_multi_nccl(int myRank,splitPoints *splitPoint, double power,int max_points_on_device,int *globalToLocalIndex,int **globalToGhostIndex,transferPoints **receiveBuffer,int *partVector,transferPoints **sendBuffer)
{
    int  k, r, nbh;
    double x_i, y_i, z_i, x_k, y_k, z_k;
    double delx, dely, delz, dist, weights;
    double sum_delx_sqr, sum_dely_sqr, sum_delz_sqr;
    double sum_delx_dely, sum_dely_delz, sum_delz_delx;
    double sum_delx_delq[5] = {0}, sum_dely_delq[5] = {0}, sum_delz_delq[5] = {0};
    double det, one_by_det;
    double temp[5], qtilde_i[5], qtilde_nbh[5];

    int bx = blockIdx.x;
    int tx = threadIdx.x;
    int i = bx * blockDim.x + tx;
    if (i < 0 || i >= max_points_on_device)
    {
        return;
    }

    x_i = splitPoint[i].x;
    y_i = splitPoint[i].y;
    z_i = splitPoint[i].z;
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
    for (k = 0; k < splitPoint[i].numberOfLocalNbhs; k++)
    //
    {
        int globalIndex = splitPoint[i].localNbhs[k];
        nbh=globalToLocalIndex[globalIndex];
        //
        x_k = splitPoint[nbh].x;
        y_k = splitPoint[nbh].y;
        z_k = splitPoint[nbh].z;
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
            qtilde_i[r] = splitPoint[i].q[r] - 0.50 * (delx * splitPoint[i].dq[0][r] + dely * splitPoint[i].dq[1][r] + delz * splitPoint[i].dq[2][r]);
            qtilde_nbh[r] = splitPoint[nbh].q[r] - 0.50 * (delx * splitPoint[nbh].dq[0][r] + dely * splitPoint[nbh].dq[1][r]+ delz * splitPoint[nbh].dq[2][r]);
            temp[r] = qtilde_nbh[r] - qtilde_i[r];
            // if(splitPoint[i].globalIndex==19595){
            // printf("Local nbh %d\n",globalIndex);
            // printf("temp[%d] is :%.15f \n",r,temp[r]);
            // printf("qtilde_i[%d] is :%.15f \n",r,qtilde_i[r]);
            // printf("qtilde_nbh[%d] is :%.15f \n",r,qtilde_nbh[r]);
            // }
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
    for (k = 0; k < splitPoint[i].numberOfGhostNbhs; k++)
    //
    {
        nbh = splitPoint[i].ghostNbhs[k];
        int device=partVector[nbh];
        int ghostIndex=globalToGhostIndex[device][nbh];
        //
        x_k = receiveBuffer[device][ghostIndex].x;
        y_k = receiveBuffer[device][ghostIndex].y;
        z_k = receiveBuffer[device][ghostIndex].z;
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
            qtilde_i[r] = splitPoint[i].q[r] - 0.50 * (delx * splitPoint[i].dq[0][r] + dely * splitPoint[i].dq[1][r] + delz * splitPoint[i].dq[2][r]);
            qtilde_nbh[r] = receiveBuffer[device][ghostIndex].q[r] - 0.50 * (delx * receiveBuffer[device][ghostIndex].dq[0][r] + dely * receiveBuffer[device][ghostIndex].dq[1][r] + delz * receiveBuffer[device][ghostIndex].dq[2][r]);
            temp[r] = qtilde_nbh[r] - qtilde_i[r];
            // if(splitPoint[i].globalIndex==19595){
            // printf("Ghost nbh %d\n",nbh);
            // printf("temp[%d] is :%.15f \n",r,temp[r]);
            // printf("qtilde_i[%d] is :%.15f \n",r,qtilde_i[r]);
            // printf("qtilde_nbh[%d] is :%.15f \n",r,qtilde_nbh[r]);
            // }
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
    // if(splitPoint[i].globalIndex==19595)
    //     printf("%d det is :%.15f \n",i,det);
    one_by_det = 1.00 / det;
    //
    for (r = 0; r < 5; r++)
    {
        temp[r] = sum_delx_delq[r] * (sum_dely_sqr * sum_delz_sqr - sum_dely_delz * sum_dely_delz) - sum_dely_delq[r] * (sum_delx_dely * sum_delz_sqr - sum_delz_delx * sum_dely_delz) + sum_delz_delq[r] * (sum_delx_dely * sum_dely_delz - sum_delz_delx * sum_dely_sqr);
    }
    //
    for (r = 0; r < 5; r++)
    {
        splitPoint[i].temp[0][r] = temp[r] * one_by_det;
        //
    }
    for (r = 0; r < 5; r++)
    {
        temp[r] = sum_delx_sqr * (sum_dely_delq[r] * sum_delz_sqr - sum_dely_delz * sum_delz_delq[r]) - sum_delx_dely * (sum_delx_delq[r] * sum_delz_sqr - sum_delz_delx * sum_delz_delq[r]) + sum_delz_delx * (sum_delx_delq[r] * sum_dely_delz - sum_delz_delx * sum_dely_delq[r]);
    }
    //
    for (r = 0; r < 5; r++)
    {
        splitPoint[i].temp[1][r] = temp[r] * one_by_det;
    }
    //
    for (r = 0; r < 5; r++)
    {
        temp[r] = sum_delx_sqr * (sum_dely_sqr * sum_delz_delq[r] - sum_dely_delq[r] * sum_dely_delz) - sum_delx_dely * (sum_delx_dely * sum_delz_delq[r] - sum_delx_delq[r] * sum_dely_delz) + sum_delz_delx * (sum_delx_dely * sum_dely_delq[r] - sum_delx_delq[r] * sum_dely_sqr);
    }
    //
    for (r = 0; r < 5; r++)
    {
        splitPoint[i].temp[2][r] = temp[r] * one_by_det;
    }
    //
}

__global__ void update_inner_loop_multi_nccl(int myRank,splitPoints *splitPoint,int max_points_on_device,transferPoints **sendBuffer)
{
    int bx = blockIdx.x;
    int tx = threadIdx.x;
    int i = bx * blockDim.x + tx;
    int r;
    if (i < 0 || i >= max_points_on_device)
    {
        return;
    }
    for (r = 0; r < 5; r++)
    {
        splitPoint[i].dq[0][r] = splitPoint[i].temp[0][r];
        splitPoint[i].dq[1][r] = splitPoint[i].temp[1][r];
        splitPoint[i].dq[2][r] = splitPoint[i].temp[2][r];
    }
    // if(i==5150 && myRank==1){
    //     for(int r=0;r<5;r++){
    //         printf("q[%d] is :%.15f \n",i,splitPoint[i].q[r]);
    //         printf("dq[0][%d] is :%.15f \n",i,splitPoint[i].dq[0][r]);
    //         printf("dq[1][%d] is :%.15f \n",i,splitPoint[i].dq[1][r]);
    //         printf("dq[2][%d] is :%.15f \n",i,splitPoint[i].dq[2][r]);
    //     }
    // }
    if(splitPoint[i].isGhost){
		for(int t=0;t<splitPoint[i].numberOfPartitionsToSendTo;++t){
			for(int r=0;r<5;++r){
				sendBuffer[splitPoint[i].partitions[t]][splitPoint[i].ghostIndex[t]].dq[0][r]=splitPoint[i].dq[0][r];
                sendBuffer[splitPoint[i].partitions[t]][splitPoint[i].ghostIndex[t]].dq[1][r]=splitPoint[i].dq[1][r];
				sendBuffer[splitPoint[i].partitions[t]][splitPoint[i].ghostIndex[t]].dq[2][r]=splitPoint[i].dq[2][r];
			}
		}
	}
}