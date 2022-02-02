int beadIndex(int index, int npc)
{

	int j,k;

	j = (int)(index/npc);

	k = j*npc;

	return (index - k);
}
