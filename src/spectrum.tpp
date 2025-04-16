// spectrum.tpp

namespace Spectrum {

template <typename Op>
inline PowerSpec scalar_arith(const PowerSpec& spec, const double scalar, Op op) {
    if (spec.is_scalar()) {
        return PowerSpec(spec.k(), op(spec.P(), scalar));
    } else {
        std::vector<double> result(spec.kvec().size());
        std::transform(spec.Pvec().begin(), spec.Pvec().end(), result.begin(), 
            [&](double p) { return op(p, scalar); });
        return PowerSpec(spec.kvec(), result);
    }
}

template <typename Op>
inline PowerSpec spec_arith(const PowerSpec& spec1, const PowerSpec& spec2, Op op) {
    if (spec1.is_scalar() && spec2.is_scalar()) {
        if (spec1.k() != spec2.k())
            throw std::invalid_argument("PowerSpec: Mismatched k values for scalar objects!");
        return PowerSpec(spec1.k(), op(spec1.P(), spec2.P()));
    } 
    else if (!spec1.is_scalar() && !spec2.is_scalar()) {
        const auto& kvec1 = spec1.kvec();
        const auto& kvec2 = spec2.kvec();
        const auto& Pvec1 = spec1.Pvec();
        const auto& Pvec2 = spec2.Pvec();

        if (kvec1.size() != kvec2.size())
            throw std::invalid_argument("PowerSpec: Mismatched k-vector sizes!");

        for (size_t i = 0; i < kvec1.size(); ++i) {
            if (kvec1[i] != kvec2[i]) {
                throw std::invalid_argument("PowerSpec: Mismatch in k-vector values at index " + std::to_string(i) + "!");
            }
        }

        std::vector<double> Pvec_new;
        Pvec_new.reserve(Pvec1.size());

        for (size_t i = 0; i < Pvec1.size(); ++i)
            Pvec_new.push_back(op(Pvec1[i], Pvec2[i]));

        return PowerSpec(kvec1, Pvec_new);
    } 
    else {
        throw std::invalid_argument("PowerSpec: Cannot mix scalar and vector types!");
    }
}

} // namespace Spectrum