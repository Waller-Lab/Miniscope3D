function y = tv3d_iso_Haar(x, tau, alpha)

% Private functions here
% circshift does circular shifting
% indexing: x(5:10), 1 indexed. Use x(5:4:end-6) to index in strides of 4
% to the 6th-to-last element
    D = 3;
    gamma = 1;   %step size
    thresh = sqrt(2) * 2 * D * tau * gamma;
    y = zeros(size(x), 'like', x);
    for axis = 1 : 3
        if axis == 3
            t_scale = alpha;
        else
            t_scale = 1;
        end
        y = y + iht3(ht3(x, axis, false, thresh*t_scale), axis, false);
        y = y + iht3(ht3(x, axis, true, thresh*t_scale), axis, true);
    end
    y = y / (2 * D);
return 

function w = ht3(x, ax, shift, thresh)
    s = size(x);
    w = zeros(s, 'like', x);
    C = 1 / sqrt(2);
    if shift
        x = circshift(x, -1, ax);
    end
    m = floor(s(ax) / 2);
    if ax == 1
        w(1:m, :, :) = C * (x(2:2:end, :, :) + x(1:2:end, :, :));  % use diff or circhisft?
        w((m + 1):end, :, :) = hs_soft(C * (x(2:2:end, :, :) - x(1:2:end, :, :)), thresh);
        %w((m + 1):end, :, :) = hs_soft(w((m + 1):end, :, :), thresh);
    elseif ax == 2
        w(:, 1:m, :) = C * (x(:, 2:2:end, :) + x(:, 1:2:end, :));
        w(:, (m + 1):end, :) = hs_soft(C * (x(:, 2:2:end, :) - x(:, 1:2:end, :)), thresh);
        %w(:, (m + 1):end, :) = hs_soft(w(:, (m + 1):end, :), thresh);
    else
        w(:, :, 1:m) = C * (x(:, :, 2:2:end) + x(:, :, 1:2:end));
        w(:, :, (m + 1):end) = hs_soft(C * (x(:, :, 2:2:end) - x(:, :, 1:2:end)), thresh);
        %w(:, :, (m + 1):end) = hs_soft(w(:, :, (m + 1):end), thresh);
    end
return

function y = iht3(w, ax, shift)
    s = size(w);
    y = zeros(s, 'like', w);
    C = 1 / sqrt(2);
    m = floor(s(ax) / 2);
    if ax == 1
        y(1:2:end, :, :) = C * (w(1:m, :, :) - w((m + 1):end, :, :));
        y(2:2:end, :, :) = C * (w(1:m, :, :) + w((m + 1):end, :, :));
    elseif ax == 2
        y(:, 1:2:end, :) = C * (w(:, 1:m, :) - w(:, (m + 1):end, :));
        y(:, 2:2:end, :) = C * (w(:, 1:m, :) + w(:, (m + 1):end, :));
    else
        y(:, :, 1:2:end) = C * (w(:, :, 1:m) - w(:, :, (m + 1):end));
        y(:, :, 2:2:end) = C * (w(:, :, 1:m) + w(:, :, (m + 1):end));
    end
    if shift
        y = circshift(y, 1, ax);
    end
return

function threshed = hs_soft(x,tau)

    threshed = max(abs(x)-tau,0);
    threshed = threshed.*sign(x);
return