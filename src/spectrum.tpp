// spectrum.tpp

namespace Spectrum {

template <typename Op>
inline PowerSpec scalar_arith(const PowerSpec& spec, const double scalar, Op op) {
    std::vector<double> result(spec.k().size());
    std::transform(spec.P().begin(), spec.P().end(), result.begin(), 
        [&](double p) { return op(p, scalar); });
    return PowerSpec(spec.k(), result);
}

template <typename Op>
inline PowerSpec spec_arith(const PowerSpec& spec1, const PowerSpec& spec2, Op op) {
    const auto& kvec1 = spec1.k();
    const auto& kvec2 = spec2.k();
    const auto& Pvec1 = spec1.P();
    const auto& Pvec2 = spec2.P();

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

} // namespace Spectrum