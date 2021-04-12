
# Convert to degrees (if necessary)
degrees(x::UN.Quantity) = convert(typeof(one(x.val)°), x)
degrees(x) = error("Please specify angles in radians or degrees. See documentation for detail")

# Convert to radians (if necessary)
radians(x::UN.Quantity) = convert(typeof(one(x.val)rad), x)
radians(x) = error("Please specify angles in radians or degrees. See documentation for detail")