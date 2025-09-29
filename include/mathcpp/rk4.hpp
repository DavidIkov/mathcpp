
#pragma once

#include <utility>

namespace mathcpp {
template <typename TimeT, typename ValueT, typename ContextT>
class RK4 {
public:
    template <typename TimeTT, typename ValueTT, typename ContextTT>
    RK4(TimeTT&& t0, ValueTT&& y0, ContextTT&& context,
        ValueT (*f)(TimeT const&, ValueT const&, ContextT const&))
        : t_(std::forward<TimeTT&&>(t0)),
          y_(std::forward<ValueTT&&>(y0)),
          context_(std::forward<ContextTT&&>(context)),
          f_(f) {}

    ValueT const& Step(TimeT const& dt) {
        ValueT k1 = f_(t_, y_, context_);
        ValueT k2 = f_(t_ + dt / 2, y_ + dt * k1 / 2, context_);
        ValueT k3 = f_(t_ + dt / 2, y_ + dt * k2 / 2, context_);
        ValueT k4 = f_(t_ + dt, y_ + dt * k3, context_);
        y_ += dt / 6 * (k1 + k2 * 2 + k3 * 2 + k4);
        t_ += dt;
        return y_;
    }

    [[nodiscard]] TimeT const& GetT() const { return t_; }
    [[nodiscard]] ValueT const& GetY() const { return y_; }
    [[nodiscard]] auto GetF() const { return f_; }

private:
    TimeT t_;
    ValueT y_;
    ContextT context_;
    ValueT (*f_)(TimeT const&, ValueT const&, ContextT const&);
};
}  // namespace mathcpp
