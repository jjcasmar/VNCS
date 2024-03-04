#ifndef VNCS_GENERATORS_BLENDINGFIELD_H
#define VNCS_GENERATORS_BLENDINGFIELD_H

namespace VNCS
{
namespace Generators
{
template <typename T, typename P>
class BlendingField
{
public:
    BlendingField() = default;
    virtual ~BlendingField() = default;

    virtual T blending(const P &p) const = 0;
};

}  // namespace Generators
}  // namespace VNCS

#endif  // VNCS_BLENDINGFIELD_H
