
#pragma once

#include <utility>

namespace mathcpp {
template <typename TimeT, typename ValueT, typename GlobalContextT,
          typename StepContextT>
class RK4 {
public:
    template <typename TimeTT, typename ValueTT, typename GlobalContextTT>
    constexpr RK4(TimeTT&& t0, ValueTT&& y0, GlobalContextTT&& global_context,
                  ValueT (*f)(TimeT const&, ValueT const&,
                              GlobalContextT const&, StepContextT const&))
        : t_(std::forward<TimeTT&&>(t0)),
          y_(std::forward<ValueTT&&>(y0)),
          global_context_(std::forward<GlobalContextTT&&>(global_context)),
          f_(f) {}

    constexpr ValueT const& Step(TimeT const& dt,
                                 StepContextT const& step_context) {
        ValueT k1 = f_(t_, y_, global_context_, step_context);
        ValueT k2 =
            f_(t_ + dt / 2, y_ + dt * k1 / 2, global_context_, step_context);
        ValueT k3 =
            f_(t_ + dt / 2, y_ + dt * k2 / 2, global_context_, step_context);
        ValueT k4 = f_(t_ + dt, y_ + dt * k3, global_context_, step_context);
        y_ += dt / 6 * (k1 + k2 * 2 + k3 * 2 + k4);
        t_ += dt;
        return y_;
    }

    [[nodiscard]] constexpr TimeT const& GetT() const { return t_; }
    [[nodiscard]] constexpr ValueT const& GetY() const { return y_; }
    [[nodiscard]] constexpr auto GetF() const { return f_; }

private:
    TimeT t_;
    ValueT y_;
    GlobalContextT global_context_;
    ValueT (*f_)(TimeT const&, ValueT const&, GlobalContextT const&,
                 StepContextT const&);
};
}  // namespace mathcpp
