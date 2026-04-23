#ifndef VIS_H
#define VIS_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(__GNUC__) || defined(__clang__)
#define VIS_INLINE static __inline__
#else
#define VIS_INLINE static
#endif

#if defined(__has_include)
#if __has_include(<X11/Xlib.h>) && __has_include(<X11/Xutil.h>)
#define VIS_HAS_X11 1
#endif
#endif

#ifndef VIS_HAS_X11
#define VIS_HAS_X11 0
#endif

enum {
    VIS_FIELD_RHO = 0,
    VIS_FIELD_SPEED = 1
};

typedef struct {
    int enabled;
    int alive;
    int lx;
    int ly;
    int scale;
#if VIS_HAS_X11
    void *display_ptr;
    unsigned long window;
    void *gc_ptr;
    void *image_ptr;
    unsigned long wm_delete;
#endif
} vis_t;

#if VIS_HAS_X11

#include <X11/Xlib.h>
#include <X11/Xutil.h>

VIS_INLINE unsigned long vis_color(double t) {
    double r, g, b;
    unsigned long ur, ug, ub;

    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;

    r = t;
    g = 1.0 - fabs(2.0 * t - 1.0);
    b = 1.0 - t;

    ur = (unsigned long)(255.0 * r);
    ug = (unsigned long)(255.0 * g);
    ub = (unsigned long)(255.0 * b);

    return (ur << 16) | (ug << 8) | ub;
}

VIS_INLINE int vis_poll(vis_t *vis) {
    Display *dpy;

    if (!vis || !vis->enabled || !vis->alive) return 0;

    dpy = (Display *)vis->display_ptr;
    while (XPending(dpy)) {
        XEvent ev;
        XNextEvent(dpy, &ev);
        if (ev.type == ClientMessage &&
            (Atom)ev.xclient.data.l[0] == (Atom)vis->wm_delete) {
            vis->alive = 0;
        }
        if (ev.type == DestroyNotify || ev.type == KeyPress) {
            vis->alive = 0;
        }
    }

    return vis->alive;
}

VIS_INLINE int vis_init(vis_t *vis, int lx, int ly, int scale, const char *title) {
    Display *dpy;
    int screen;
    int width, height;
    XImage *img;

    if (!vis) return 0;
    memset(vis, 0, sizeof(*vis));

    if (lx <= 0 || ly <= 0) return 0;
    if (scale < 1) scale = 1;

    dpy = XOpenDisplay(NULL);
    if (!dpy) {
        fprintf(stderr, "vis: cannot open X11 display\n");
        return 0;
    }

    screen = DefaultScreen(dpy);
    width = lx * scale;
    height = ly * scale;

    vis->window = XCreateSimpleWindow(
        dpy, RootWindow(dpy, screen), 0, 0,
        (unsigned int)width, (unsigned int)height, 1,
        BlackPixel(dpy, screen), WhitePixel(dpy, screen));
    XStoreName(dpy, (Window)vis->window, title ? title : "LBM");
    XSelectInput(dpy, (Window)vis->window, ExposureMask | KeyPressMask | StructureNotifyMask);

    vis->wm_delete = (unsigned long)XInternAtom(dpy, "WM_DELETE_WINDOW", False);
    XSetWMProtocols(dpy, (Window)vis->window, (Atom *)&vis->wm_delete, 1);
    XMapWindow(dpy, (Window)vis->window);

    img = XCreateImage(
        dpy, DefaultVisual(dpy, screen), (unsigned int)DefaultDepth(dpy, screen),
        ZPixmap, 0, NULL, width, height, 32, 0);
    if (!img) {
        XDestroyWindow(dpy, (Window)vis->window);
        XCloseDisplay(dpy);
        fprintf(stderr, "vis: cannot create X11 image\n");
        return 0;
    }

    img->data = (char *)calloc((size_t)img->bytes_per_line * (size_t)img->height, 1);
    if (!img->data) {
        img->data = NULL;
        XDestroyImage(img);
        XDestroyWindow(dpy, (Window)vis->window);
        XCloseDisplay(dpy);
        fprintf(stderr, "vis: out of memory\n");
        return 0;
    }

    vis->enabled = 1;
    vis->alive = 1;
    vis->lx = lx;
    vis->ly = ly;
    vis->scale = scale;
    vis->display_ptr = dpy;
    vis->gc_ptr = (void *)DefaultGC(dpy, screen);
    vis->image_ptr = img;

    XFlush(dpy);
    return 1;
}

VIS_INLINE double vis_scalar(const double *rholb, const double *uxlb, const double *uylb, int lx, int ly, int f, int x, int y, int field) {
    double rho;

    rho = rholb[x+lx*(ly-1-y)];

    if (field == VIS_FIELD_RHO) {
        return rho;
    }

    if (rho <= 1e-14) {
        return 0.0;
    }
// enum {_C, _N, _S, _W, _E, _NE, _NW, _SW, _SE};
    {
        double ux = uxlb[x+lx*(ly-1-y)];
        double uy = uylb[x+lx*(ly-1-y)];
        return sqrt(ux * ux + uy * uy);
    }
}

VIS_INLINE int vis_draw(vis_t *vis, const double *rholb, const double *uxlb, const double *uylb, int lx, int ly, int f, int field) {
    Display *dpy;
    XImage *img;
    int x, y, dx, dy;
    double vmin, vmax;

    if (!vis || !vis->enabled || !vis->alive || !rholb || !uxlb || !uylb) return 0;
    if (lx != vis->lx || ly != vis->ly || f < 5) return 0;
    if (!vis_poll(vis)) return 0;

    vmin = vis_scalar(rholb, uxlb, uylb, lx, ly, f, 0, 0, field);
    vmax = vmin;
    for (y = 0; y < ly; ++y) {
        for (x = 0; x < lx; ++x) {
            double v = vis_scalar(rholb, uxlb, uylb, lx, ly, f, x, y, field);
            if (v < vmin) vmin = v;
            if (v > vmax) vmax = v;
        }
    }

    dpy = (Display *)vis->display_ptr;
    img = (XImage *)vis->image_ptr;

    for (y = 0; y < ly; ++y) {
        for (x = 0; x < lx; ++x) {
            double v = vis_scalar(rholb, uxlb, uylb, lx, ly, f, x, y, field);
            double t = (vmax > vmin) ? (v - vmin) / (vmax - vmin) : 0.0;
            unsigned long pixel = vis_color(t);

            for (dy = 0; dy < vis->scale; ++dy) {
                for (dx = 0; dx < vis->scale; ++dx) {
                    XPutPixel(img, x * vis->scale + dx, y * vis->scale + dy, pixel);
                }
            }
        }
    }

    XPutImage(dpy, (Window)vis->window, (GC)vis->gc_ptr, img, 0, 0, 0, 0,
              (unsigned int)(lx * vis->scale), (unsigned int)(ly * vis->scale));
    XFlush(dpy);

    return vis_poll(vis);
}

VIS_INLINE int vis_draw_rho(vis_t *vis, const double *rholb, const double *uxlb, const double *uylb, int lx, int ly, int f) {
    return vis_draw(vis, rholb, uxlb, uylb, lx, ly, f, VIS_FIELD_RHO);
}

VIS_INLINE int vis_draw_speed(vis_t *vis, const double *rholb, const double *uxlb, const double *uylb, int lx, int ly, int f) {
    return vis_draw(vis, rholb, uxlb, uylb, lx, ly, f, VIS_FIELD_SPEED);
}

VIS_INLINE void vis_close(vis_t *vis) {
    Display *dpy;

    if (!vis || !vis->enabled) return;

    dpy = (Display *)vis->display_ptr;
    if (vis->image_ptr) {
        XDestroyImage((XImage *)vis->image_ptr);
        vis->image_ptr = NULL;
    }
    if (vis->window) {
        XDestroyWindow(dpy, (Window)vis->window);
        vis->window = 0;
    }
    if (dpy) {
        XCloseDisplay(dpy);
    }

    memset(vis, 0, sizeof(*vis));
}

#else

VIS_INLINE int vis_poll(vis_t *vis) {
    (void)vis;
    return 0;
}

VIS_INLINE int vis_init(vis_t *vis, int lx, int ly, int scale, const char *title) {
    (void)lx;
    (void)ly;
    (void)scale;
    (void)title;
    if (vis) memset(vis, 0, sizeof(*vis));
    fprintf(stderr, "vis: X11 headers not available at compile time\n");
    return 0;
}

VIS_INLINE double vis_scalar(const double *rholb, const double *uxlb, const double *uylb, int lx, int ly, int f, int x, int y, int field) {
    (void)cells;
    (void)lx;
    (void)ly;
    (void)f;
    (void)x;
    (void)y;
    (void)field;
    return 0.0;
}

VIS_INLINE int vis_draw(vis_t *vis, const double *rholb, const double *uxlb, const double *uylb, int lx, int ly, int f, int field) {
    (void)vis;
    (void)rholb;
    (void)uxlb;
    (void)uylb;
    (void)lx;
    (void)ly;
    (void)f;
    (void)field;
    return 0;
}

VIS_INLINE int vis_draw_rho(vis_t *vis, const double *rholb, const double *uxlb, const double *uylb, int lx, int ly, int f) {
    (void)vis;
    (void)rholb;
    (void)uxlb;
    (void)uylb;
    (void)lx;
    (void)ly;
    (void)f;
    return 0;
}

VIS_INLINE int vis_draw_speed(vis_t *vis, const double *rholb, const double *uxlb, const double *uylb, int lx, int ly, int f) {
    (void)vis;
    (void)rholb;
    (void)uxlb;
    (void)uylb;
    (void)lx;
    (void)ly;
    (void)f;
    return 0;
}

    (void)lx;
    (void)ly;
    (void)f;
    return 0;
}

VIS_INLINE int vis_draw_speed(vis_t *vis, const double *rholb, const double *uxlb, const double *uylb, int lx, int ly, int f) {
    (void)vis;
    (void)rholb;
    (void)uxlb;
    (void)uylb;
    (void)lx;
    (void)ly;
    (void)f;
    return 0;
}

VIS_INLINE void vis_close(vis_t *vis) {
    (void)vis;
}

#endif

#endif
