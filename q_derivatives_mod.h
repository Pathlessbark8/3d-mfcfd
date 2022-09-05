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
            point.qm[i][0][k] = point.q[i][k]; // q_maximum ..
            point.qm[i][1][k] = point.q[i][k]; // q_minimum ..
        }
        for (k = 0; k < point.nbhs[i]; k++)
        //
        {
            nbh = point.conn[i][k];
            //
            for (r = 0; r < 5; r++)
            {
                if (point.q[nbh][r] > point.qm[i][0][r])
                {
                    point.qm[i][0][r] = point.q[nbh][r];
                }
                if (point.q[nbh][r] < point.qm[i][1][r])
                {
                    point.qm[i][1][r] = point.q[nbh][r];
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
                temp[r] = (point.q[nbh][r] - point.q[i][r]);

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
            point.dq[i][0][k] = temp[k] * one_by_det;
        }
        //
        for (k = 0; k < 5; k++)
        {
            temp[k] = sum_delx_sqr * (sum_dely_delq[k] * sum_delz_sqr - sum_dely_delz * sum_delz_delq[k]) - sum_delx_dely * (sum_delx_delq[k] * sum_delz_sqr - sum_delz_delx * sum_delz_delq[k]) + sum_delz_delx * (sum_delx_delq[k] * sum_dely_delz - sum_delz_delx * sum_dely_delq[k]);
        }
        //
        for (k = 0; k < 5; k++)
        {
            point.dq[i][1][k] = temp[k] * one_by_det;
        }
        //
        for (k = 0; k < 5; k++)
        {
            temp[k] = sum_delx_sqr * (sum_dely_sqr * sum_delz_delq[k] - sum_dely_delq[k] * sum_dely_delz) - sum_delx_dely * (sum_delx_dely * sum_delz_delq[k] - sum_delx_delq[k] * sum_dely_delz) + sum_delz_delx * (sum_delx_dely * sum_dely_delq[k] - sum_delx_delq[k] * sum_dely_sqr);
        }
        for (k = 0; k < 5; k++)
        {
            point.dq[i][2][k] = temp[k] * one_by_det;
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
            nbh = point.conn[i][k];
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
                qtilde_i[r] = point.q[i][r] - 0.50 * (delx * point.dq[i][0][r] + dely * point.dq[i][1][r] + delz * point.dq[i][2][r]);
                qtilde_nbh[r] = point.q[nbh][r] - 0.50 * (delx * point.dq[nbh][0][r] + dely * point.dq[nbh][1][r] + delz * point.dq[nbh][2][r]);
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
            point.temp[i][0][r] = temp[r] * one_by_det;
            //
        }
        for (r = 0; r < 5; r++)
        {
            temp[r] = sum_delx_sqr * (sum_dely_delq[r] * sum_delz_sqr - sum_dely_delz * sum_delz_delq[r]) - sum_delx_dely * (sum_delx_delq[r] * sum_delz_sqr - sum_delz_delx * sum_delz_delq[r]) + sum_delz_delx * (sum_delx_delq[r] * sum_dely_delz - sum_delz_delx * sum_dely_delq[r]);
        }
        //
        for (r = 0; r < 5; r++)
        {
            point.temp[i][1][r] = temp[r] * one_by_det;
        }
        //
        for (r = 0; r < 5; r++)
        {
            temp[r] = sum_delx_sqr * (sum_dely_sqr * sum_delz_delq[r] - sum_dely_delq[r] * sum_dely_delz) - sum_delx_dely * (sum_delx_dely * sum_delz_delq[r] - sum_delx_delq[r] * sum_dely_delz) + sum_delz_delx * (sum_delx_dely * sum_dely_delq[r] - sum_delx_delq[r] * sum_dely_sqr);
        }
        //
        for (r = 0; r < 5; r++)
        {
            point.temp[i][2][r] = temp[r] * one_by_det;
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
                point.dq[i][0][r] = point.temp[i][0][r];
                point.dq[i][1][r] = point.temp[i][1][r];
                point.dq[i][2][r] = point.temp[i][2][r];
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
        point.qm[i][0][k] = point.q[i][k]; // q_maximum ..
        point.qm[i][1][k] = point.q[i][k]; // q_minimum ..
    }
    for (k = 0; k < point.nbhs[i]; k++)
    //
    {
        nbh = point.conn[i][k];
        //
        for (r = 0; r < 5; r++)
        {
            if (point.q[nbh][r] > point.qm[i][0][r])
            {
                point.qm[i][0][r] = point.q[nbh][r];
            }
            if (point.q[nbh][r] < point.qm[i][1][r])
            {
                point.qm[i][1][r] = point.q[nbh][r];
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
            temp[r] = (point.q[nbh][r] - point.q[i][r]);

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
        point.dq[i][0][k] = temp[k] * one_by_det;
    }
    //
    for (k = 0; k < 5; k++)
    {
        temp[k] = sum_delx_sqr * (sum_dely_delq[k] * sum_delz_sqr - sum_dely_delz * sum_delz_delq[k]) - sum_delx_dely * (sum_delx_delq[k] * sum_delz_sqr - sum_delz_delx * sum_delz_delq[k]) + sum_delz_delx * (sum_delx_delq[k] * sum_dely_delz - sum_delz_delx * sum_dely_delq[k]);
    }
    //
    for (k = 0; k < 5; k++)
    {
        point.dq[i][1][k] = temp[k] * one_by_det;
    }
    //
    for (k = 0; k < 5; k++)
    {
        temp[k] = sum_delx_sqr * (sum_dely_sqr * sum_delz_delq[k] - sum_dely_delq[k] * sum_dely_delz) - sum_delx_dely * (sum_delx_dely * sum_delz_delq[k] - sum_delx_delq[k] * sum_dely_delz) + sum_delz_delx * (sum_delx_dely * sum_dely_delq[k] - sum_delx_delq[k] * sum_dely_sqr);
    }
    for (k = 0; k < 5; k++)
    {
        point.dq[i][2][k] = temp[k] * one_by_det;
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
        nbh = point.conn[i][k];
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
            qtilde_i[r] = point.q[i][r] - 0.50 * (delx * point.dq[i][0][r] + dely * point.dq[i][1][r] + delz * point.dq[i][2][r]);
            qtilde_nbh[r] = point.q[nbh][r] - 0.50 * (delx * point.dq[nbh][0][r] + dely * point.dq[nbh][1][r] + delz * point.dq[nbh][2][r]);
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
        point.temp[i][0][r] = temp[r] * one_by_det;
        //
    }
    for (r = 0; r < 5; r++)
    {
        temp[r] = sum_delx_sqr * (sum_dely_delq[r] * sum_delz_sqr - sum_dely_delz * sum_delz_delq[r]) - sum_delx_dely * (sum_delx_delq[r] * sum_delz_sqr - sum_delz_delx * sum_delz_delq[r]) + sum_delz_delx * (sum_delx_delq[r] * sum_dely_delz - sum_delz_delx * sum_dely_delq[r]);
    }
    //
    for (r = 0; r < 5; r++)
    {
        point.temp[i][1][r] = temp[r] * one_by_det;
    }
    //
    for (r = 0; r < 5; r++)
    {
        temp[r] = sum_delx_sqr * (sum_dely_sqr * sum_delz_delq[r] - sum_dely_delq[r] * sum_dely_delz) - sum_delx_dely * (sum_delx_dely * sum_delz_delq[r] - sum_delx_delq[r] * sum_dely_delz) + sum_delz_delx * (sum_delx_dely * sum_dely_delq[r] - sum_delx_delq[r] * sum_dely_sqr);
    }
    //
    for (r = 0; r < 5; r++)
    {
        point.temp[i][2][r] = temp[r] * one_by_det;
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
        point.dq[i][0][r] = point.temp[i][0][r];
        point.dq[i][1][r] = point.temp[i][1][r];
        point.dq[i][2][r] = point.temp[i][2][r];
    }
}
